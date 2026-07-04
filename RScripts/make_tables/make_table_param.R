# ==============================================================================
# Prepare table showing model parameterization
# ==============================================================================

m <- Model.final$m
q <- Model.final$q

decimal <- 3

#c@{\\extracolsep{\\fill}}rrrrrrr

latex_table <- rbind("\\begin{table}[H]",
                     "\\caption{Model parameterization}",
                     "\\label{tab:param}",
                     "\\renewcommand{\\arraystretch}{0.7}",
                     "\\begin{tabular*}{\\textwidth}",
                     "{r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}r@{\\extracolsep{\\fill}}}",
                     "\\hline",
                     "\\multicolumn{8}{c}{{\\bf Panel A - Trend and cycle loadings}}\\\\",
                     "\\hline",
                     "&$\\delta^{(\\pi)}_C$ & $\\delta^{(\\Delta y)}_C$ & & $\\delta^{(\\pi)}_T$ & $\\delta^{(\\Delta y)}_T$\\\\",
                     "\\hline")

for(i in 1:m){
  this_line <- paste("$\\delta_{C_", i, "}$ &",
                     make.entry(Model.final$delta.c[i,1],decimal),"&",
                     make.entry(Model.final$delta.c[i,2],decimal),"&",
                     "$\\delta_{T_", i, "}$ &",
                     make.entry(Model.final$delta.t[i,1],decimal),"&",
                     make.entry(Model.final$delta.t[i,2],decimal),"\\\\",
                     sep="")
  latex_table <- rbind(latex_table,
                       this_line)
}
latex_table <- rbind(latex_table,
                     "\\hline")

latex_table <- rbind(latex_table,
                     "\\multicolumn{8}{c}{{\\bf Panel B - Dynamics of $\\mathcal{Y}_t$}}\\\\",
                     "\\hline",
                     "\\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\Phi_\\mathcal{Y}$} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\Theta$} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\Gamma_{\\mathcal{Y},0}$} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\Gamma_{\\mathcal{Y},1}$}\\\\",
                     "\\hline")

line1 <- paste("$\\phi_{1,1}$ &",
               make.entry(Model.final$Phi.Y.r[1,1],decimal),
               "& $\\Theta^s$ &",
               make.entry(Model.final$Theta[1,1],decimal),
               "& $\\Gamma_{1,\\mathcal{Y},0}$ &",
               make.entry(Model.final$Gamma.Y0.r[1,1],decimal),
               "& $\\Gamma_{[3,5],\\mathcal{Y},1}$ &",
               make.entry(Model.final$Gamma.Y1.r[5,3],decimal),"\\\\",
               sep="")

line2 <- paste("$\\phi_{2,2}$ &",
               make.entry(Model.final$Phi.Y.r[2,2],decimal),
               "&$\\Theta^d$ &",
               make.entry(Model.final$Theta[2,3],decimal),
               "&$\\Gamma_{2,\\mathcal{Y},0}$ &",
               make.entry(Model.final$Gamma.Y0.r[2,1],decimal),
               "& $\\Gamma_{[4,5],\\mathcal{Y},1}$ &",
               make.entry(Model.final$Gamma.Y1.r[5,4],decimal),"\\\\",
               sep="")
line3 <- paste("$\\phi_{3,3}$ &",
               make.entry(Model.final$Phi.Y.r[3,3],decimal),
               "& &",
               "&$\\Gamma_{3,\\mathcal{Y},0}$ &",
               make.entry(Model.final$Gamma.Y0.r[3,1],decimal),
               "& &","\\\\",
               sep="")
line4 <- paste("$\\phi_{4,4}$ &",
               make.entry(Model.final$Phi.Y.r[4,4],decimal),
               "& &",
               "&$\\Gamma_{4,\\mathcal{Y},0}$ &",
               make.entry(Model.final$Gamma.Y0.r[4,1],decimal),
               "& &","\\\\",
               "\\hline",
               sep="")

latex_table <- rbind(latex_table,
                     line1, line2, line3, line4
)


latex_table <- rbind(latex_table,
                     "\\multicolumn{8}{c}{{\\bf Panel C - Dynamics of $z_t$}}\\\\",
                     "\\hline",
                     "\\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\nu$} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\phi$} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{$\\mu$} \\\\",
                     "\\hline")

for(i in 1:q){
  this_line <- paste("$\\nu_{", i, "}$ &",
                     make.entry(Model.final$nu[i,1],decimal),
                     "& $\\phi_{", i, i,"}$ &",
                     make.entry(Model.final$phi[i,i],decimal),
                     "& $\\mu_{", i, "}$ &",
                     make.entry(Model.final$mu[i,1],decimal),"\\\\",
                     sep="")
  latex_table <- rbind(latex_table,
                       this_line)
}



latex_table <- rbind(latex_table,
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     "\\begin{spacing}{0.8}",
                     paste("\\parbox{\\linewidth}{\\textit{Notes}: This table shows the parameter estimates. We also have: $\\rho^{\\pi}= ",
                           make.entry(Model.final$pi.bar[1],decimal,dollar=0),"$",
                           " and $\\rho^{\\Delta y}= ",
                           make.entry(Model.final$pi.bar[2],decimal,dollar=0),"$. The first panel shows the trend and cycle loadings. The second shows the dynamics of $\\mathcal{Y}_t$. The last panel shows the dynamics of $z_t$}",
                           sep = ""),
                     "\\end{spacing}",
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_param"
latex.file <- paste(name.of.file,".txt", sep="")
if (!exists("path_table")) {
  if (exists("area") && area == "EA") {
    path_table <- "tables/EA_2024/Baseline/"
  } else {
    path_table <- "tables/US_2024/Baseline/"
  }
}
dir.create(path_table, recursive = TRUE, showWarnings = FALSE)
write(latex_table, file.path(path_table, latex.file))
