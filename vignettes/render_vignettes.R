  # utility for making vignette
  # given an Rmd markdown filename prefix from the vignettes directory,
  # create the md and html files
  
  # you could create an array of vignette names and loop the rendering
  Rmd_name <- "STIRexample"
  # Rmd render function (below) is run relative to the STIR/ director
  input1 <- paste("vignettes/",Rmd_name,".Rmd",sep="")      
  # render once to keep (no clean) md file, and save the rendering object (temp)
  temp <- rmarkdown::render(input1, run_pandoc = FALSE, clean = FALSE)
  knit_meta <- attr(temp, "knit_meta") 
  # render again creates html if you want that
  Rmd_name <- "STIRexample"  # running the Rmd seems to lose Rmd_name variable
  input2 <- paste("vignettes/",Rmd_name,".knit.md",sep="")
  rmarkdown::render(input = input2, knit_meta = knit_meta )
  
  # when you render this way, it creates variables in the console environment
  rm(list=ls())
