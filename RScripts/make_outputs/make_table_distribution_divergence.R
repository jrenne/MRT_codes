# ==============================================================================
# Prepare table showing distribution divergence
# ==============================================================================

decimal <- 3

horizon <- 1

### 1YEAR
#for(horizon in c(1)){
latex_table <- rbind("\\begin{table}[ph!]",
                     "\\begingroup \\tiny",
                     paste("\\caption{Distribution Divergence Comparison Table}"),
                     "\\label{tab:distribution.divergence}",
                     paste("\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rrrrrr}",sep=""),
                     "\\hline",
                     "\\hline",
                     " \\multicolumn{2}{c}{a) Inflation} \\\\",
                     " & \\multicolumn{2}{c}{Baseline} & \\multicolumn{2}{c}{No Kurtosis} & \\multicolumn{2}{c}{No Higher Order Moments} \\\\",
                     "\\cmidrule(lr){2-3}", "\\cmidrule(lr){4-5}", "\\cmidrule(lr){6-7}",
                     " & $d_{TV}$ & $d_{KL}$ & $d_{TV}$ & $d_{KL}$ & $d_{TV}$ & $d_{KL}$ \\\\",
                     "\\hline"
                     )



#### INFLATION
for(horizon in 5:8){
  
latex_table <- rbind(latex_table,
                     paste(" \\multicolumn{1}{c}{horizon = ", horizon,"Q} \\\\"),
                     "\\cmidrule(lr){1-1}"
)

eval(parse(text = gsub(" ","",paste("TVD.4th <- TVD.",horizon,"Q.4th", sep=""))))
eval(parse(text = gsub(" ","",paste("KLD.4th <- KLD.",horizon,"Q.4th", sep=""))))
eval(parse(text = gsub(" ","",paste("WD.4th <- WD.",horizon,"Q.4th", sep=""))))

eval(parse(text = gsub(" ","",paste("TVD.with.3rd.only <- TVD.",horizon,"Q.with.3rd.only", sep=""))))
eval(parse(text = gsub(" ","",paste("KLD.with.3rd.only <- KLD.",horizon,"Q.with.3rd.only", sep=""))))
eval(parse(text = gsub(" ","",paste("WD.with.3rd.only <- WD.",horizon,"Q.with.3rd.only", sep=""))))

eval(parse(text = gsub(" ","",paste("TVD.no.3rd.4th <- TVD.",horizon,"Q.no.3rd.4th", sep=""))))
eval(parse(text = gsub(" ","",paste("KLD.no.3rd.4th <- KLD.",horizon,"Q.no.3rd.4th", sep=""))))
eval(parse(text = gsub(" ","",paste("WD.no.3rd.4th <- WD.",horizon,"Q.no.3rd.4th", sep=""))))


this_line <- paste("Min, Max","&(",
                   make.entry(min(TVD.4th[,2]),decimal),", ", make.entry(max(TVD.4th[,2]),decimal),")&(",
                   make.entry(min(KLD.4th[,2]),decimal),", ", make.entry(max(KLD.4th[,2]),decimal),")&(",
                   #make.entry(min(WD.4th[,2]),decimal),", ", make.entry(max(WD.4th[,2]),decimal),")&(",
                   
                   make.entry(min(TVD.with.3rd.only[,2]),decimal),", ", make.entry(max(TVD.with.3rd.only[,2]),decimal),")&(",
                   make.entry(min(KLD.with.3rd.only[,2]),decimal),", ", make.entry(max(KLD.with.3rd.only[,2]),decimal),")&(",
                   #make.entry(min(WD.with.3rd.only[,2]),decimal),", ", make.entry(max(WD.with.3rd.only[,2]),decimal),")&(",
                   
                   make.entry(min(TVD.no.3rd.4th[,2]),decimal),", ", make.entry(max(TVD.no.3rd.4th[,2]),decimal),")&(",
                   make.entry(min(KLD.no.3rd.4th[,2]),decimal),", ", make.entry(max(KLD.no.3rd.4th[,2]),decimal),")",
                   #make.entry(min(WD.no.3rd.4th[,2]),decimal),", ", make.entry(max(WD.no.3rd.4th[,2]),decimal),")",
                   "\\\\",
                   sep="")

this_line2 <- paste("Median","&",
                   make.entry(median(TVD.4th[,2]),decimal),"&",
                   make.entry(median(KLD.4th[,2]),decimal),"&",
                   #make.entry(median(WD.4th[,2]),decimal),"&",
                   
                   make.entry(median(TVD.with.3rd.only[,2]),decimal),"&",
                   make.entry(median(KLD.with.3rd.only[,2]),decimal),"&",
                   #make.entry(median(WD.with.3rd.only[,2]),decimal),"&",
                   
                   make.entry(median(TVD.no.3rd.4th[,2]),decimal),"&",
                   make.entry(median(KLD.no.3rd.4th[,2]),decimal),"&",
                   #make.entry(median(WD.no.3rd.4th[,2]),decimal),
                   "\\\\",
                   sep="")


this_line1 <- paste("25th, 75th","&(",
                    make.entry(quantile(TVD.4th[,2],0.25),decimal),", ", make.entry(quantile(TVD.4th[,2],0.75),decimal),")&(",
                    make.entry(quantile(KLD.4th[,2],0.25),decimal),", ", make.entry(quantile(KLD.4th[,2],0.75),decimal),")&(",
                    #make.entry(quantile(WD.4th[,2],0.25),decimal),", ", make.entry(quantile(WD.4th[,2],0.75),decimal),")&(",
                    
                    make.entry(quantile(TVD.with.3rd.only[,2],0.25),decimal),", ", make.entry(quantile(TVD.with.3rd.only[,2],0.75),decimal),")&(",
                    make.entry(quantile(KLD.with.3rd.only[,2],0.25),decimal),", ", make.entry(quantile(KLD.with.3rd.only[,2],0.75),decimal),")&(",
                    #make.entry(quantile(WD.with.3rd.only[,2],0.25),decimal),", ", make.entry(quantile(WD.with.3rd.only[,2],0.75),decimal),")&(",
                    
                    make.entry(quantile(TVD.no.3rd.4th[,2],0.25),decimal),", ", make.entry(quantile(TVD.no.3rd.4th[,2],0.75),decimal),")&(",
                    make.entry(quantile(KLD.no.3rd.4th[,2],0.25),decimal),", ", make.entry(quantile(KLD.no.3rd.4th[,2],0.75),decimal),")",
                    #make.entry(quantile(WD.no.3rd.4th[,2],0.25),decimal),", ", make.entry(quantile(WD.no.3rd.4th[,2],0.75),decimal),")",
                    "\\\\",
                    sep="")


latex_table <- rbind(latex_table,
                     this_line,
                     this_line1,
                     this_line2,
                     "\\\\"
)


}

latex_table <- rbind(latex_table,
                     paste(" \\multicolumn{1}{c}{Overall", "","} \\\\"),
                     "\\cmidrule(lr){1-1}"
)

TVD.4th <- c(TVD.5Q.4th[,2],TVD.6Q.4th[,2],TVD.7Q.4th[,2],TVD.8Q.4th[,2])
KLD.4th <- c(KLD.5Q.4th[,2],KLD.6Q.4th[,2],KLD.7Q.4th[,2],KLD.8Q.4th[,2])
WD.4th <- c(WD.5Q.4th[,2],WD.6Q.4th[,2],WD.7Q.4th[,2],WD.8Q.4th[,2])

TVD.with.3rd.only <- c(TVD.5Q.with.3rd.only[,2],TVD.6Q.with.3rd.only[,2],TVD.7Q.with.3rd.only[,2],TVD.8Q.with.3rd.only[,2])
KLD.with.3rd.only <- c(KLD.5Q.with.3rd.only[,2],KLD.6Q.with.3rd.only[,2],KLD.7Q.with.3rd.only[,2],KLD.8Q.with.3rd.only[,2])
WD.with.3rd.only <- c(WD.5Q.with.3rd.only[,2],WD.6Q.with.3rd.only[,2],WD.7Q.with.3rd.only[,2],KLD.8Q.with.3rd.only[,2])

TVD.no.3rd.4th <- c(TVD.5Q.no.3rd.4th[,2],TVD.6Q.no.3rd.4th[,2],TVD.7Q.no.3rd.4th[,2],TVD.8Q.no.3rd.4th[,2])
KLD.no.3rd.4th <- c(KLD.5Q.no.3rd.4th[,2],KLD.6Q.no.3rd.4th[,2],KLD.7Q.no.3rd.4th[,2],KLD.8Q.no.3rd.4th[,2])
WD.no.3rd.4th <- c(WD.5Q.no.3rd.4th[,2],WD.6Q.no.3rd.4th[,2],WD.7Q.no.3rd.4th[,2],KLD.8Q.no.3rd.4th[,2])

this_line <- paste("Min, Max","&(",
                   make.entry(min(TVD.4th),decimal),", ", make.entry(max(TVD.4th),decimal),")&(",
                   make.entry(min(KLD.4th),decimal),", ", make.entry(max(KLD.4th),decimal),")&(",
                   #make.entry(min(WD.4th),decimal),", ", make.entry(max(WD.4th),decimal),")&(",
                   
                   make.entry(min(TVD.with.3rd.only),decimal),", ", make.entry(max(TVD.with.3rd.only),decimal),")&(",
                   make.entry(min(KLD.with.3rd.only),decimal),", ", make.entry(max(KLD.with.3rd.only),decimal),")&(",
                   #make.entry(min(WD.with.3rd.only),decimal),", ", make.entry(max(WD.with.3rd.only),decimal),")&(",
                   
                   make.entry(min(TVD.no.3rd.4th),decimal),", ", make.entry(max(TVD.no.3rd.4th),decimal),")&(",
                   make.entry(min(KLD.no.3rd.4th),decimal),", ", make.entry(max(KLD.no.3rd.4th),decimal),")",
                   #make.entry(min(KLD.no.3rd.4th),decimal),", ", make.entry(max(KLD.no.3rd.4th),decimal),")",
                   "\\\\",
                   sep="")

this_line2 <- paste("Median","&",
                    make.entry(median(TVD.4th),decimal),"&",
                    make.entry(median(KLD.4th),decimal),"&",
                    #make.entry(median(WD.4th),decimal),"&",
                    
                    make.entry(median(TVD.with.3rd.only),decimal),"&",
                    make.entry(median(KLD.with.3rd.only),decimal),"&",
                    #make.entry(median(WD.with.3rd.only),decimal),"&",
                    
                    make.entry(median(TVD.no.3rd.4th),decimal),"&",
                    make.entry(median(KLD.no.3rd.4th),decimal),"",
                    #make.entry(median(WD.no.3rd.4th),decimal),
                    "\\\\",
                    sep="")

this_line1 <- paste("25th, 75th","&(",
                    make.entry(quantile(TVD.4th,0.25),decimal),", ", make.entry(quantile(TVD.4th,0.75),decimal),")&(",
                    make.entry(quantile(KLD.4th,0.25),decimal),", ", make.entry(quantile(KLD.4th,0.75),decimal),")&(",
                    #make.entry(quantile(WD.4th,0.25),decimal),", ", make.entry(quantile(WD.4th,0.75),decimal),")&(",
                    
                    make.entry(quantile(TVD.with.3rd.only,0.25),decimal),", ", make.entry(quantile(TVD.with.3rd.only,0.75),decimal),")&(",
                    make.entry(quantile(KLD.with.3rd.only,0.25),decimal),", ", make.entry(quantile(KLD.with.3rd.only,0.75),decimal),")&(",
                    #make.entry(quantile(WD.with.3rd.only,0.25),decimal),", ", make.entry(quantile(WD.with.3rd.only,0.75),decimal),")&(",
                    
                    make.entry(quantile(TVD.no.3rd.4th,0.25),decimal),", ", make.entry(quantile(TVD.no.3rd.4th,0.75),decimal),")&(",
                    make.entry(quantile(KLD.no.3rd.4th,0.25),decimal),", ", make.entry(quantile(KLD.no.3rd.4th,0.75),decimal),")",
                    #make.entry(quantile(WD.no.3rd.4th,0.25),decimal),", ", make.entry(quantile(KLD.no.3rd.4th,0.75),decimal),")",
                    "\\\\",
                    sep="")


latex_table <- rbind(latex_table,
                     this_line,
                     this_line1,
                     this_line2,
                     "\\\\"
)


latex_table <- rbind(latex_table,
                     "\\hline",
                     " \\multicolumn{2}{c}{b) Real GDP growth} \\\\",
                     "\\\\"
                     )


#### GDP
for(horizon in 5:8){
  
  latex_table <- rbind(latex_table,
                       paste(" \\multicolumn{1}{c}{horizon = ", horizon,"Q} \\\\"),
                       "\\cmidrule(lr){1-1}"
  )
  
  eval(parse(text = gsub(" ","",paste("TVD.4th <- TVD.",horizon,"Q.4th.gdp", sep=""))))
  eval(parse(text = gsub(" ","",paste("KLD.4th <- KLD.",horizon,"Q.4th.gdp", sep=""))))
  eval(parse(text = gsub(" ","",paste("WD.4th <- WD.",horizon,"Q.4th.gdp", sep=""))))
  
  eval(parse(text = gsub(" ","",paste("TVD.with.3rd.only <- TVD.",horizon,"Q.with.3rd.only.gdp", sep=""))))
  eval(parse(text = gsub(" ","",paste("KLD.with.3rd.only <- KLD.",horizon,"Q.with.3rd.only.gdp", sep=""))))
  eval(parse(text = gsub(" ","",paste("WD.with.3rd.only <- WD.",horizon,"Q.with.3rd.only.gdp", sep=""))))
  
  eval(parse(text = gsub(" ","",paste("TVD.no.3rd.4th <- TVD.",horizon,"Q.no.3rd.4th.gdp", sep=""))))
  eval(parse(text = gsub(" ","",paste("KLD.no.3rd.4th <- KLD.",horizon,"Q.no.3rd.4th.gdp", sep=""))))
  eval(parse(text = gsub(" ","",paste("WD.no.3rd.4th <- WD.",horizon,"Q.no.3rd.4th.gdp", sep=""))))
  
  
  this_line <- paste("Min, Max","&(",
                     make.entry(min(TVD.4th[,2]),decimal),", ", make.entry(max(TVD.4th[,2]),decimal),")&(",
                     make.entry(min(KLD.4th[,2]),decimal),", ", make.entry(max(KLD.4th[,2]),decimal),")&(",
                     #make.entry(min(WD.4th[,2]),decimal),", ", make.entry(max(WD.4th[,2]),decimal),")&(",
                     
                     make.entry(min(TVD.with.3rd.only[,2]),decimal),", ", make.entry(max(TVD.with.3rd.only[,2]),decimal),")&(",
                     make.entry(min(KLD.with.3rd.only[,2]),decimal),", ", make.entry(max(KLD.with.3rd.only[,2]),decimal),")&(",
                     #make.entry(min(WD.with.3rd.only[,2]),decimal),", ", make.entry(max(WD.with.3rd.only[,2]),decimal),")&(",
                     
                     make.entry(min(TVD.no.3rd.4th[,2]),decimal),", ", make.entry(max(TVD.no.3rd.4th[,2]),decimal),")&(",
                     make.entry(min(KLD.no.3rd.4th[,2]),decimal),", ", make.entry(max(KLD.no.3rd.4th[,2]),decimal),")",
                     #make.entry(min(WD.no.3rd.4th[,2]),decimal),", ", make.entry(max(WD.no.3rd.4th[,2]),decimal),")",
                     
                     "\\\\",
                     sep="")
  
  this_line2 <- paste("Median","&",
                      make.entry(median(TVD.4th[,2]),decimal),"&",
                      make.entry(median(KLD.4th[,2]),decimal),"&",
                      #make.entry(median(WD.4th[,2]),decimal),"&",
                      
                      make.entry(median(TVD.with.3rd.only[,2]),decimal),"&",
                      make.entry(median(KLD.with.3rd.only[,2]),decimal),"&",
                      #make.entry(median(WD.with.3rd.only[,2]),decimal),"&",
                      
                      make.entry(median(TVD.no.3rd.4th[,2]),decimal),"&",
                      make.entry(median(KLD.no.3rd.4th[,2]),decimal),"",
                      #make.entry(median(WD.no.3rd.4th[,2]),decimal),
                      "\\\\",
                      sep="")
  
  this_line1 <- paste("25th, 75th","&(",
                     make.entry(quantile(TVD.4th[,2],0.25),decimal),", ", make.entry(quantile(TVD.4th[,2],0.75),decimal),")&(",
                     make.entry(quantile(KLD.4th[,2],0.25),decimal),", ", make.entry(quantile(KLD.4th[,2],0.75),decimal),")&(",
                     #make.entry(quantile(WD.4th[,2],0.25),decimal),", ", make.entry(quantile(WD.4th[,2],0.75),decimal),")&(",
                     
                     make.entry(quantile(TVD.with.3rd.only[,2],0.25),decimal),", ", make.entry(quantile(TVD.with.3rd.only[,2],0.75),decimal),")&(",
                     make.entry(quantile(KLD.with.3rd.only[,2],0.25),decimal),", ", make.entry(quantile(KLD.with.3rd.only[,2],0.75),decimal),")&(",
                     #make.entry(quantile(WD.with.3rd.only[,2],0.25),decimal),", ", make.entry(quantile(WD.with.3rd.only[,2],0.75),decimal),")&(",
                     
                     make.entry(quantile(TVD.no.3rd.4th[,2],0.25),decimal),", ", make.entry(quantile(TVD.no.3rd.4th[,2],0.75),decimal),")&(",
                     make.entry(quantile(KLD.no.3rd.4th[,2],0.25),decimal),", ", make.entry(quantile(KLD.no.3rd.4th[,2],0.75),decimal),")",
                     #make.entry(quantile(WD.no.3rd.4th[,2],0.25),decimal),", ", make.entry(quantile(WD.no.3rd.4th[,2],0.75),decimal),")",
                     "\\\\",
                     sep="")
  
  
  latex_table <- rbind(latex_table,
                       this_line,
                       this_line1,
                       this_line2,
                       "\\\\"
  )
  
}

latex_table <- rbind(latex_table,
                     paste(" \\multicolumn{1}{c}{Overall", "","} \\\\"),
                     "\\cmidrule(lr){1-1}"
)

TVD.4th <- c(TVD.5Q.4th.gdp[,2],TVD.6Q.4th.gdp[,2],TVD.7Q.4th.gdp[,2],TVD.8Q.4th.gdp[,2])
KLD.4th <- c(KLD.5Q.4th.gdp[,2],KLD.6Q.4th.gdp[,2],KLD.7Q.4th.gdp[,2],KLD.8Q.4th.gdp[,2])
WD.4th <- c(WD.5Q.4th.gdp[,2],WD.6Q.4th.gdp[,2],WD.7Q.4th.gdp[,2],WD.8Q.4th.gdp[,2])

TVD.with.3rd.only <- c(TVD.5Q.with.3rd.only.gdp[,2],TVD.6Q.with.3rd.only.gdp[,2],TVD.7Q.with.3rd.only.gdp[,2],TVD.8Q.with.3rd.only.gdp[,2])
KLD.with.3rd.only <- c(KLD.5Q.with.3rd.only.gdp[,2],KLD.6Q.with.3rd.only.gdp[,2],KLD.7Q.with.3rd.only.gdp[,2],KLD.8Q.with.3rd.only.gdp[,2])
WD.with.3rd.only <- c(WD.5Q.with.3rd.only.gdp[,2],WD.6Q.with.3rd.only.gdp[,2],WD.7Q.with.3rd.only.gdp[,2],WD.8Q.with.3rd.only.gdp[,2])

TVD.no.3rd.4th <- c(TVD.5Q.no.3rd.4th.gdp[,2],TVD.6Q.no.3rd.4th.gdp[,2],TVD.7Q.no.3rd.4th.gdp[,2],TVD.8Q.no.3rd.4th.gdp[,2])
KLD.no.3rd.4th <- c(KLD.5Q.no.3rd.4th.gdp[,2],KLD.6Q.no.3rd.4th.gdp[,2],KLD.7Q.no.3rd.4th.gdp[,2],KLD.8Q.no.3rd.4th.gdp[,2])
WD.no.3rd.4th <- c(WD.5Q.no.3rd.4th.gdp[,2],WD.6Q.no.3rd.4th.gdp[,2],WD.7Q.no.3rd.4th.gdp[,2],WD.8Q.no.3rd.4th.gdp[,2])

this_line <- paste("Min, Max","&(",
                   make.entry(min(TVD.4th),decimal),", ", make.entry(max(TVD.4th),decimal),")&(",
                   make.entry(min(KLD.4th),decimal),", ", make.entry(max(KLD.4th),decimal),")&(",
                   #make.entry(min(WD.4th),decimal),", ", make.entry(max(WD.4th),decimal),")&(",
                   
                   make.entry(min(TVD.with.3rd.only),decimal),", ", make.entry(max(TVD.with.3rd.only),decimal),")&(",
                   make.entry(min(KLD.with.3rd.only),decimal),", ", make.entry(max(KLD.with.3rd.only),decimal),")&(",
                   #make.entry(min(WD.with.3rd.only),decimal),", ", make.entry(max(WD.with.3rd.only),decimal),")&(",
                   
                   make.entry(min(TVD.no.3rd.4th),decimal),", ", make.entry(max(TVD.no.3rd.4th),decimal),")&(",
                   make.entry(min(KLD.no.3rd.4th),decimal),", ", make.entry(max(KLD.no.3rd.4th),decimal),")",
                   #make.entry(min(WD.no.3rd.4th),decimal),", ", make.entry(max(WD.no.3rd.4th),decimal),")",
                   
                   "\\\\",
                   sep="")

this_line2 <- paste("Median","&",
                    make.entry(median(TVD.4th),decimal),"&",
                    make.entry(median(KLD.4th),decimal),"&",
                    #make.entry(median(WD.4th),decimal),"&",
                    
                    make.entry(median(TVD.with.3rd.only),decimal),"&",
                    make.entry(median(KLD.with.3rd.only),decimal),"&",
                    #make.entry(median(WD.with.3rd.only),decimal),"&",
                    
                    make.entry(median(TVD.no.3rd.4th),decimal),"&",
                    make.entry(median(KLD.no.3rd.4th),decimal),"",
                    #make.entry(median(WD.no.3rd.4th),decimal),
                    "\\\\",
                    sep="")

this_line1 <- paste("25th, 75th","&(",
                    make.entry(quantile(TVD.4th,0.25),decimal),", ", make.entry(quantile(TVD.4th,0.75),decimal),")&(",
                    make.entry(quantile(KLD.4th,0.25),decimal),", ", make.entry(quantile(KLD.4th,0.75),decimal),")&(",
                    #make.entry(quantile(WD.4th,0.25),decimal),", ", make.entry(quantile(WD.4th,0.75),decimal),")&(",
                    
                    make.entry(quantile(TVD.with.3rd.only,0.25),decimal),", ", make.entry(quantile(TVD.with.3rd.only,0.75),decimal),")&(",
                    make.entry(quantile(KLD.with.3rd.only,0.25),decimal),", ", make.entry(quantile(KLD.with.3rd.only,0.75),decimal),")&(",
                    #make.entry(quantile(WD.with.3rd.only,0.25),decimal),", ", make.entry(quantile(WD.with.3rd.only,0.75),decimal),")&(",
                    
                    make.entry(quantile(TVD.no.3rd.4th,0.25),decimal),", ", make.entry(quantile(TVD.no.3rd.4th,0.75),decimal),")&(",
                    make.entry(quantile(KLD.no.3rd.4th,0.25),decimal),", ", make.entry(quantile(KLD.no.3rd.4th,0.75),decimal),")",
                    #make.entry(quantile(WD.no.3rd.4th,0.25),decimal),", ", make.entry(quantile(WD.no.3rd.4th,0.75),decimal),")",
                    
                    "\\\\",
                    sep="")


latex_table <- rbind(latex_table,
                     this_line,
                     this_line1,
                     this_line2,
                     "\\\\"
)


latex_table <- rbind(latex_table,
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{\\textit{Notes}:",
                     "}",
                     "\\end{footnotesize}",
                     "\\endgroup",
                     "\\end{table}")


name.of.file <- "table_distribution_divergence"
latex.file <- paste(name.of.file,".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))

