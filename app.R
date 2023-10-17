#1st function mz vs RKMD####
library(rsconnect)
library(shiny)
library(tidyverse)
library("xlsx")
library(ggplot2)
library(bslib)
library(colourpicker) #new
library(gridExtra) #new
library(ggpubr) #new
arrangeKMD_Rcode <- function(datafile, a, b, c, d, e, f, g, h5) {
  datafile -> tibbleA
  a->a_init
  NA -> tibbleA[8][tibbleA[9] < d]
  NA -> tibbleA[8][tibbleA[9] > e]
  #Mass Intensity Filter, where d = lower limit; e = upper limit
  #IMPORTANT! NA's is NOT 0; NA's will NOT be filtered out
  
  tibbleA[[1]][complete.cases(tibbleA[[1]])] -> RKMDName_noNA
  tibbleA[[2]][complete.cases(tibbleA[[2]])] -> RKMDTheoMass_noNA
  tibbleA[[8]][complete.cases(tibbleA[[8]])] -> MassList_noNA
  #pulls column out of tibble A and returns as a vector without NA's
  #g is repeating unit exact mass (14.01565 for CH2)
  #h is the saturation mass (e.g. H, F, Cl, Br)
  
  RKMDCalc <- function(x, y) {
    if ((floor(x) %% 2) == (floor(y) %% 2)) {
      (((((x * (floor(g) / g)) - (floor(x * (floor(g) / g)))
      ) - ((y * (floor(g) / g))
           - (floor(y * (floor(g) /  g)))
      ))) / (h5*2*(floor(g) / g) - 2*floor(h5)))
    }
    else {
      print(NA)
    }
  }
  #function for calculating REFERENCED KMD, x is exp. mass list and y is RKMD theo. mass
  #will also filter out even/odd masses if it doesn't match the theo Reference mass
  
  
  tibble(title = 1:length(RKMDTheoMass_noNA)) -> tibbleB
  for (i in seq_along(MassList_noNA)) {
    unlist(MassList_noNA[i] %>% map2(RKMDTheoMass_noNA, RKMDCalc)) -> tibbleB[i]
  }
  #map2 to calculate new lists from two different vectors
  #arguments that vary come BEFORE function
  
  paste(MassList_noNA) -> MassListString
  MassListString -> colnames(tibbleB)
  
  #filtering parameters: a,b,c
  #a-b is range of interger, a MUST be lower than b
  #c is % to filter, e.g. c=10%, then keep 0.9-1.1, 1.9-2.1, etc.
  
  NA -> tibbleB[tibbleB <= (a - (c / 100))]
  NA -> tibbleB[tibbleB >= (b + (c / 100))]
  
  while (a < b) {
    NA -> tibbleB[tibbleB >= (a + (c / 100)) & tibbleB <= (a + 1 - (c / 100))]
    {
      a = a + 1
    }
  }
  a_init->a
#Radyl carbon filter
  RC_Calc <- function(x, y, z) {  
    abs(x-(y-2.0157*z))/g
  }
  
  tibble(title = 1:length(RKMDTheoMass_noNA)) -> tibble_RC_Filter
  
  for (i in seq_along(MassList_noNA)) {
    unlist(pmap(list(RKMDTheoMass_noNA,MassList_noNA[i],tibbleB[[i]]), RC_Calc)) -> tibble_RC_Filter[i]
  }
  
  -100->x
  100->y
  if (f=="Both"){
    while (x < y) {
      NA -> tibbleB[tibble_RC_Filter >= (x + (c / 100)) & tibble_RC_Filter <= (x + 1 - (c / 100))]
      {
        x = x + 1
      }
    } 
  }else if (f=="Even"){
    
    while (x < y) {
      NA -> tibbleB[tibble_RC_Filter >= (x + (c / 100)) & tibble_RC_Filter <= (x + 1 - (c / 100)) 
                    | round(tibble_RC_Filter)%%2==1]
      {
        x = x + 1
      }
    }
  } else if (f=="Odd"){
    
    while (x < y) {
      NA -> tibbleB[tibble_RC_Filter >= (x + (c / 100)) & tibble_RC_Filter <= (x + 1 - (c / 100))
                    | round(tibble_RC_Filter)%%2==0]
      {
        x = x + 1
      }
    }
  } else if (f=="None") { } 
  
  tibbleB %>% round(3)%>% mutate(RKMDNames = RKMDName_noNA, .before = 1) -> tibbleB
  tibbleB[,colSums(is.na(tibbleB))<nrow(tibbleB)]-> tibbleB 
  #remove all columns with all NA
  #tibbleB is a calculated data frame of RKMD, NA's applied to conditions
  #all empty columns removed
  
  rep(list(1), ncol(tibbleB)) -> listA
  for (i in seq_along(tibbleB)) {
    tibbleB %>%
      arrange(tibbleB[[i]]) -> listA[[i]]
    listA[[i]] %>%
      filter(!is.na(listA[[i]][i])) -> listA[[i]]
  }
  
  nrow(tibbleB) -> h
  1 -> i
  
  while (i <= length(listA)) {
    if (nrow(listA[[i]]) == h) {
      i + 1 -> i
    } else {
      add_row(listA[[i]]) -> listA[[i]]
    }
  }
  print(listA)
  #arranges each column by i , makes a list of tibbles (one for each arranged column)
  #also filters out NA's in i column
  #Adds empty rows to make all tibbles in list have the same dimensions
  
  
  rep(list(1), length(listA)) -> listB
  tibble(1:(nrow(listA[[1]]))) -> df1
  for (i in seq_along(listA)) {
    pull(listA[[i]], var = 1) %>% tibble() -> df1[1]
    pull(listA[[i]], var = i) %>% tibble() -> df1[2]
    df1 -> listB[[i]]
  }
  print(listB)
  #listB:pulls out the 1st and i column to match the relevant arranged columns
  
  tibble(1:nrow(tibbleB)) -> outputA
  for (i in seq_along(1:length(listB))) {
    listB[[i]][[1]] -> outputA[[(i * 2) - 1]]
    listB[[i]][[2]] -> outputA[[(i * 2)]]
  }
  
  #outputA is listB with the items in the list going across instead of down
  
  colnames(tibbleB) -> names
  tibble(1:nrow(tibbleB), 1) -> vectorname
  
  c("name1") -> vectorname
  for (i in seq_along(1:ncol(outputA))) {
    if (i %% 2 == 0) {
      names[(i / 2)] -> vectorname[i]
    } else {
      "RKMDName" -> vectorname[i]
    }
  }
  vectorname -> colnames(outputA)
  #rename the columns correctly
  return(outputA)
  
}
#2nd function#######
arrangeKMD_RKMDvsmz_count_repeats <- function(datafile, a, b, c, d, e, f, g, h5) {
  a -> varA
  datafile -> tibbleA
  
  NA -> tibbleA[8][tibbleA[9] < d]
  NA -> tibbleA[8][tibbleA[9] > e]
  #Mass Intensity Filter, where d = lower limit; e = upper limit
  #IMPORTANT! NA's is NOT 0; NA's will NOT be filtered out
  
  tibbleA[[1]][complete.cases(tibbleA[[1]])] -> RKMDName_noNA
  tibbleA[[2]][complete.cases(tibbleA[[2]])] -> RKMDTheoMass_noNA
  tibbleA[[8]][complete.cases(tibbleA[[8]])] -> MassList_noNA
  #pulls column out of tibble A and returns as a vector without NA's
  
  
  RKMDCalc <- function(x, y) {
    if ((floor(x) %% 2) == (floor(y) %% 2)) {
      (((((x * (floor(g) / g)) - (floor(x * (floor(g) / g)))
      ) - ((y * (floor(g) / g))
           - (floor(y * (floor(g) /  g)))
      ))) /(h5*2*(floor(g) / g) - 2*floor(h5)))
    }
    else {
      print(NA)
    }
  } 
  #function for calculating REFERENCED KMD, x is exp. mass list and y is RKMD theo. mass
  #will also filter out even/odd masses if it doesn't match the theo Reference mass
  
  
  tibble(title = 1:length(RKMDTheoMass_noNA)) -> tibbleB
  for (i in seq_along(MassList_noNA)) {
    unlist(MassList_noNA[i] %>% map2(RKMDTheoMass_noNA, RKMDCalc)) -> tibbleB[i]
  }
  #map2 to calculate new lists from two different vectors
  #arguments that vary come BEFORE function
  
  paste(MassList_noNA) -> MassListString
  MassListString -> colnames(tibbleB)
  

  #filtering parameters: a,b,c
  #a-b is range of interger, a MUST be lower than b
  #c is % to filter, e.g. c=10%, then keep 0.9-1.1, 1.9-2.1, etc.
  
  NA -> tibbleB[tibbleB <= (a - (c / 100))]
  NA -> tibbleB[tibbleB >= (b + (c / 100))]
  while (a < b) {
    NA -> tibbleB[tibbleB >= (a + (c / 100)) & tibbleB <= (a + 1 - (c / 100))]
    {
      a = a + 1
    }
  }
  varA->a
  RC_Calc <- function(x, y, z) {  
  abs(x-(y-2.0157*z))/g
  }
  
  tibble(title = 1:length(RKMDTheoMass_noNA)) -> tibble_RC_Filter
  
  for (i in seq_along(MassList_noNA)) {
    unlist(pmap(list(RKMDTheoMass_noNA,MassList_noNA[i],tibbleB[[i]]), RC_Calc)) -> tibble_RC_Filter[i]
  }
  
-100->x
100->y
if (f=="Both"){
  while (x < y) {
    NA -> tibbleB[tibble_RC_Filter >= (x + (c / 100)) & tibble_RC_Filter <= (x + 1 - (c / 100))]
    {
    x = x + 1
    }
  } 
}else if (f=="Even"){

  while (x < y) {
  NA -> tibbleB[(tibble_RC_Filter >= (x + (c / 100)) & tibble_RC_Filter <= (x + 1 - (c / 100)))
                | round(tibble_RC_Filter)%%2==1]
    {
    x = x + 1
    }
  }
} else if (f=="Odd"){

while (x < y) {
  NA -> tibbleB[(tibble_RC_Filter >= (x + (c / 100)) & tibble_RC_Filter <= (x + 1 - (c / 100)))
               | round(tibble_RC_Filter)%%2==0]
    {
    x = x + 1
    }
  }
} else if (f=="None") { } 
  
  tibbleB %>% round(3)%>% mutate(RKMDNames = RKMDName_noNA, .before = 1) -> tibbleB
  tibbleB[,colSums(is.na(tibbleB))<nrow(tibbleB)]-> tibbleB 
  #remove all columns with all NA
  #tibbleB is a calculated data frame of RKMD, NA's applied to conditions
  #all empty columns removed
  
  rep(list(1), ncol(tibbleB)) -> listA
  for (i in seq_along(tibbleB)) {
    tibbleB %>%
      arrange(tibbleB[[i]]) -> listA[[i]]
    listA[[i]] %>%
      filter(!is.na(listA[[i]][i])) -> listA[[i]]
  }
  
  nrow(tibbleB) -> h
  1 -> i
  
  while (i <= length(listA)) {
    if (nrow(listA[[i]]) == h) {
      i + 1 -> i
    } else {
      add_row(listA[[i]]) -> listA[[i]]
    }
  }
  #arranges each column by i , makes a list of tibbles (one for each arranged column)
  #also filters out NA's in i column
  #Adds empty rows to make all tibbles in list have the same dimensions
  
  rep(list(1), length(listA)) -> listB
  tibble(1:(nrow(listA[[1]]))) -> df1
  for (i in seq_along(listA)) {
    pull(listA[[i]], var = 1) %>% tibble() -> df1[1]
    pull(listA[[i]], var = i) %>% tibble() -> df1[2]
    df1 -> listB[[i]]
  }
  #listB:pulls out the 1st and i column to match the relevant arranged columns
  
  tibble(1:nrow(tibbleB)) -> outputA
  for (i in seq_along(1:length(listB))) {
    listB[[i]][[1]] -> outputA[[(i * 2) - 1]]
    listB[[i]][[2]] -> outputA[[(i * 2)]]
  }
  
  #outputA is listB with the items in the list going across instead of down
  
  colnames(tibbleB) -> names
  tibble(1:nrow(tibbleB), 1) -> vectorname
  
  c("name1") -> vectorname
  for (i in seq_along(1:ncol(outputA))) {
    if (i %% 2 == 0) {
      names[(i / 2)] -> vectorname[i]
    } else {
      "RKMDName" -> vectorname[i]
    }
  }
  vectorname -> colnames(outputA)
  #rename the columns correctly
  
  tibble(1, 2) -> tibbleC
  for (i in (seq_along(1:(ncol(outputA) / 2)))) {
    length((outputA[[2 * i]])[complete.cases(outputA[[2 * i]])]) -> tibbleC[i, 2]
    as.double(colnames(outputA[2 * i])) -> tibbleC[i, 1]
  }
  tibbleC %>% filter(!is.na(tibbleC[[1]])) -> tibbleC
  
  
  #tibbleC, table of m/z values with their repeats
  
  tibble(title = 1:length(MassList_noNA)) -> tibbleB2
  for (i in seq_along(RKMDTheoMass_noNA)) {
   (-unlist(RKMDTheoMass_noNA[i] %>% map2(MassList_noNA, RKMDCalc))) -> tibbleB2[i]
  }
  tibble(title = 1:length(MassList_noNA)) -> tibble_RC_Filter2
  
  for (i in seq_along(RKMDTheoMass_noNA)) {
    unlist(pmap(list(RKMDTheoMass_noNA[i],MassList_noNA,tibbleB2[[i]]), RC_Calc)) -> tibble_RC_Filter2[i]
  }
  
  -100->x
  100->y
  if (f=="Both"){
    while (x < y) {
      NA -> tibbleB2[tibble_RC_Filter2 >= (x + (c / 100)) & tibble_RC_Filter2 <= (x + 1 - (c / 100))]
      {
        x = x + 1
      }
    } 
  }else if (f=="Even"){
    
    while (x < y) {
      NA -> tibbleB2[(tibble_RC_Filter2 >= (x + (c / 100)) & tibble_RC_Filter2 <= (x + 1 - (c / 100)))
                    | round(tibble_RC_Filter2)%%2==1]
      {
        x = x + 1
      }
    }
  } else if (f=="Odd"){
    
    while (x < y) {
      NA -> tibbleB2[(tibble_RC_Filter2 >= (x + (c / 100)) & tibble_RC_Filter2 <= (x + 1 - (c / 100)))
                    | round(tibble_RC_Filter2)%%2==0]
      {
        x = x + 1
      }
    }
  } else if (f=="None") { } 
  

  paste(RKMDName_noNA) -> MassListString
  MassListString -> colnames(tibbleB2)
  
  #tibbleB2 is RKMD vs m/z (from function #2)
  
  varA -> a
  #resets variable a (avoids the a variable change from loop earlier)
  
  NA -> tibbleB2[tibbleB2 <= (a - (c / 100))]
  NA -> tibbleB2[tibbleB2 >= (b + (c / 100))]
  while (a < b) {
    NA -> tibbleB2[tibbleB2 >= (a + (c / 100)) &
                     tibbleB2 <= (a + 1 - (c / 100))]
    {
      a = a + 1
    }
  }
  
  tibbleB2%>%round(3)->tibbleB2
  tibbleB2 %>% discard( ~ all(is.na(.))) -> tibbleB2
  tibbleB2 %>% mutate("m/z" = MassList_noNA, .before = 1) -> tibbleB2
  
  #remove all columns with all NA
  #tibbleB is a calculated data frame of RKMD, NA's applied to conditions
  
  rep(list(1), ncol(tibbleB2)) -> listA2
  for (i in seq_along(tibbleB2)) {
    tibbleB2 %>%
      arrange(tibbleB2[[i]]) -> listA2[[i]]
    listA2[[i]] %>%
      filter(!is.na(listA2[[i]][i])) -> listA2[[i]]
    listA2[[i]] %>% arrange(listA2[[i]][["m/z"]]) ->listA2[[i]]
    
  }
  
  nrow(tibbleB2) -> h
  1 -> i
  while (i <= length(listA2)) {
    if (nrow(listA2[[i]]) == h) {
      i + 1 -> i
    } else {
      add_row(listA2[[i]]) -> listA2[[i]]
    }
  }
  
  #arranges each column by i , makes a list of tibbles (one for each arranged column)
  #also filters out NA's in i column
  #Adds empty rows to make all tibbles in list have the same dimensions
  
  rep(list(1), length(listA2)) -> listrepeat
  for (i in seq_along(listA2)) {
    (tibbleC[[2]][tibbleC[[1]] %in% listA2[[i]][[1]]]) %>% tibble() -> listrepeat[[i]]
  }
  nrow(tibbleB2) -> h
  1 -> i
  while (i <= length(listrepeat)) {
    if (nrow(listrepeat[[i]]) == h) {
      i + 1 -> i
    } else {
      add_row(listrepeat[[i]]) -> listrepeat[[i]]
    }
  }
  
  #listrepeat is a list of repeats corrosponding to listA2 with the correct rows
  
  
  
  rep(list(1), length(listA2)) -> listB2
  tibble(1:(nrow(listA2[[1]]))) -> df2
  for (i in seq_along(listA2)) {
    pull(listA2[[i]], var = 1) %>% tibble() -> df2[1]
    pull(listA2[[i]], var = i) %>% tibble() -> df2[2]
    pull(listrepeat[[i]], var = 1) %>% tibble() -> df2[3]
    df2 -> listB2[[i]]
  }
  #listB2:pulls out the 1st and i column to match the relevant arranged columns
  #also adds repeats in 3rd column of each list
  
  tibble(1:nrow(tibbleB2)) -> outputA2
  for (i in seq_along(1:length(listB2))) {
    listB2[[i]][[1]] -> outputA2[[(i * 3) - 2]]
    listB2[[i]][[2]] -> outputA2[[(i * 3) - 1]]
    listB2[[i]][[3]] -> outputA2[[(i * 3)]]
    
  }
  #outputA2 is listB2 with the items in the list going across instead of down
  
  colnames(tibbleB2) -> names
  tibble(1:nrow(tibbleB2), 1) -> vectorname
  c("name1") -> vectorname
  for (i in seq_along(1:ncol(outputA2))) {
    if (i %% 3 == 2) {
      names[((i + 2) / 3)] -> vectorname[i]
    }
    if (i %% 3 == 0) {
      "repeats" -> vectorname[i]
    }
    if (i %% 3 == 1) {
      "m/z" -> vectorname[i]
    }
  }
  vectorname -> colnames(outputA2)
  #rename the columns correctly
  
  bind_cols(outputA2[2],outputA2[3]) %>% arrange("m/z") ->tibbleF1
  tibble("m/z list"=outputA2[[1]],"repeats list"=2)->tibbleF
  tibbleF1[[1,1]] %in% tibbleC[[1]]
  for(i in seq_along(outputA2[[2]])){
    if(tibbleF1[[i,1]] %in% tibbleC[[1]]) {
      ( tibbleC[,2][tibbleC[,1] == tibbleF1[[i,1]]] )->tibbleF[i,2]
    }else{0->tibbleF[i,2]}
  }
  tibbleF %>% mutate(int=datafile[[9]], .before =`repeats list`) -> tibbleF
  tibbleF[[3]] ->outputA2[3]
  return(outputA2)
  
}


#graphing RKMD plots#####

#makes a list of m/z,intensities, and repeats 
clean <- function(data) {
  c("") -> colorlist
  for (i in seq_along(1:length(data))) {
    if (data[[i]] == 1) {
      colorlist[[i]] <- "1"
    }
    else if (data[[i]] == 2) {
      colorlist[[i]] <- "2"
    }
    else if (data[[i]] == 3) {
      colorlist[[i]] <- "3"
    }
    else{
      colorlist[[i]] <- "4+"
    }
  }
  return(colorlist)
}
#Returns list of all plots
generate <-
  function(data,
           shapelist,
           colorlist,
           size,
           background1,
           lines,
           xlower,
           xupper,
           ylower,
           yupper) {
    data <- data[-(1:3)]
    list() -> plotlist
    list() -> titlelist
    
    for (i in seq_along(1:(ncol(data) / 3))) {
      m <- data[[3 * i - 2]][complete.cases(data[[3 * i - 2]])]
      md <- data[[3 * i - 1]][complete.cases(data[[3 * i - 1]])]
      c <- clean(data[[3 * i]][complete.cases(data[[3 * i]])])
      df <- data.frame(Mass = m,
                       RKMD = md,
                       Repeat = c)
      colnames(data)[[3 * i - 1]] -> titlelist[[i]]
      bracket_ind <- unlist(gregexpr("]", titlelist[[i]]))[1]
      main_str <- substr(titlelist[[i]], 0, bracket_ind)
      
      leftover <-
        substr(titlelist[[i]], bracket_ind + 2, str_length(titlelist[[i]]))
      main_str <- bquote(.(main_str) ^ "+" ~ .(leftover))
      len <- length(md)
      cols <-
        c(
          "1" = colorlist[[1]],
          "2" = colorlist[[2]],
          "3" = colorlist[[3]],
          "4+" = colorlist[[4]]
        )
      fillcols <-
        c(
          "1" = colorlist[[1]],
          "2" = colorlist[[2]],
          "3" = colorlist[[3]],
          "4+" = colorlist[[4]]
        )
      shapes <-
        c(
          "1" = shapelist[[1]],
          "2" = shapelist[[2]],
          "3" = shapelist[[3]],
          "4+" = shapelist[[4]]
        )
      
      df <-
        df %>% ggplot(aes(x = Mass, y = RKMD)) + geom_point(aes(
          colour = Repeat,
          fill = Repeat,
          shape = Repeat
        ), size = size) +
        labs(title = main_str, x = "m/z") +
        scale_shape_manual(values = shapes) +
        scale_fill_manual(values = fillcols) + scale_color_manual(values = cols) +
        #edit.. need to make this adjustable? limits and breaks 
        scale_y_continuous(limits= c(ylower,yupper), breaks = seq(ylower,yupper,by=1)) +
        scale_x_continuous(limits=c(xlower,xupper))+
                             
        #end edit
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(face = "italic"),
          panel.background = element_rect(fill = background1),
          panel.grid.major = element_line(colour = lines),
          panel.grid.minor = element_line(colour = background1)
        )
      
      #Y-axis tick marks increment by 1 unless range is too small
    
      plotlist[[titlelist[[i]]]] <- df
    }
    list(plotlist, titlelist) -> listB
    return(listB)
  }    
#plot list is the list of plots and title list is the list of titles

#graphing mz plots####
#graphs mz spectra based on repeats (color coded) based on all RKMD filters 
graphSpectrum <- function(data2,
                          colorlist2,
                          background2,
                          xaxislowerlimits2,xaxislimits2, xaxisinterval2,
                          yaxislimits2, yaxisinterval2
) {
  data2%>%ggplot()+geom_linerange(aes(x=x,y=y,ymax=y,ymin=0,color=factor(z)))+
    scale_color_manual(values = colorlist2)+
    scale_y_continuous(breaks = seq(0,yaxislimits2, by =yaxisinterval2),limits=c(0,yaxislimits2),expand = c(0, 0))+
    scale_x_continuous(breaks = seq(0,xaxislimits2, by =xaxisinterval2),limits=c(xaxislowerlimits2,xaxislimits2))+
    
    theme(panel.background=element_rect(fill=background2),
          panel.grid.major=element_line(colour=background2),
          panel.grid.minor=element_blank(),
          axis.line=element_line(colour="black"),
          axis.title.x=element_text(face="italic"))+
    
    labs(color="Legend Title")+
    xlab("m/z")+
    ylab("Intensity")
}  

#SHINY##############


#ui      #####
ui <- fluidPage(tabsetPanel(
  tabPanel(
    "Compute KMD",
    titlePanel ("KMD Analysis"),
    wellPanel(
      downloadButton("download0", "Download Template"),
      div(style = "height:20px"), #white space
      fileInput("upload", "Upload a File", accept = c(".csv", ".tsv")),
    ),
    #div(style = "height:50px"), #white space
    wellPanel(
      numericInput("inputG", "Repeating Unit Monoisotopic Mass. No Nitrogen. Default is set to CH2", value =14.01565),
      numericInput("inputH", "Repeating Unit's Saturation Monoisotopic Mass (e.g. H,F,Cl,Br).", value =1.007825),
      
      numericInput("inputA", "Lower RKMD Limit", value = -10),
      numericInput("inputB", "Upper RKMD Limit", value = 0),
      numericInput("inputC", "Percent from RKMD", value = 10),
      numericInput("intensityInputA", "Lower m/z Intensity Limit", value =
                     1),
      numericInput("intensityInputB", "Upper m/z Intensity Limit", value =
                     1000000000),
      radioButtons("RC_FilterF","Radyl Carbon Filtering for Lipids",c("None","Even","Odd","Both"))
    ),
    wellPanel(
      actionButton("button1", "Compute m/z vs RKMD; Press ONCE"),
      actionButton("button2", "Compute RKMD vs m/z; Press ONCE"),
    ),
    
    dataTableOutput("TableOutput"),
    dataTableOutput("TableOutput2"),
    
    wellPanel(
      downloadButton("download1", "Download m/z vs RKMD"),
      downloadButton("download2", "Download RKMD vs m/z")
    ),
  ),
  
  tabPanel(
    "Plot KMD",
#graphing UI    #######
    titlePanel ("KMD Plots"),
    sidebarLayout(
      sidebarPanel(
        sliderInput(
          "height",
          "Graph Height",
          min = 200,
          max = 1000,
          value = 500
        ),
        numericInput("inputxlower", "X Axis lower Limit", value =100),
        numericInput("inputxupper", "X Axis upper Limit", value =1000),
        numericInput("inputylower", "Y Axis lower Limit", value =-10),
        numericInput("inputyupper", "Y Axis upper Limit", value =1),
        sliderInput(
          "width",
          "Graph Width",
          min = 200,
          max = 1000,
          value = 500
        ),
     
        div(style = "height:22px", h4("No Overlap")),
        div(style = "display:inline-block;width:49%;height:60px;vertical-align:top;text-align:center",selectInput("shapechoices1", "Shape", c(1:25), 21)),
        div(style = "display:inline-block;width:49%;height:60px;text-align:center", colourInput("colorchoices1", "Color", value = "black")),
        div(style = "height:22px", h4("1 Overlap")),
        div(style = "display:inline-block;width:49%;height:60px;vertical-align:top;text-align:center", selectInput("shapechoices2", "Shape", c(1:25), 22)),
        div(style = "display:inline-block;width:49%;height:60px;text-align:center", colourInput("colorchoices2", "Color", value = "black")),
        div(style = "height:22px", h4("2 Overlap")),
        div(style = "display:inline-block;width:49%;height:60px;vertical-align:top;text-align:center", selectInput("shapechoices3", "Shape", c(1:25), 23)),
        div(style = "display:inline-block;width:49%;height:60px;text-align:center", colourInput("colorchoices4", "Color", value = "black")),
        div(style = "height:22px", h4("3+ Overlap")),
        div(style = "display:inline-block;width:49%;height:65px;vertical-align:top;text-align:center", selectInput("shapechoices4", "Shape", c(1:25), 24)),
        div(style = "display:inline-block;width:49%;height:65px;text-align:center", colourInput("colorchoices3", "Color", value = "black")),
        div(style = "text-align:center;height:65px", numericInput("sizechoices1", "Size of Points", value =
                                                                    2)),
        div(
          style = "display:inline-block;width:49%;text-align:center",
          colourInput("backgroundchoices1", "Background", value = "white")
        ),
        div(style = "display:inline-block;width:49%;text-align:center", colourInput("linechoices1", "Grid", value =
                                                                                      "grey")),
        plotOutput("shapeKey"),
      ),
      mainPanel(
        uiOutput("plotselect"),
        actionButton("plotButton", "Press to Plot or Apply Changes"),
        downloadButton("downloadg2", "Download RKMD vs m/z graph"),
        plotOutput(outputId = "plot1ID"),
      ),
    ),
  ),
#MS Spectra UI #######
tabPanel(
  "Mass Spectra",
  titlePanel ("Mass Spectra"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "height2",
        "Graph Height",
        min = 200,
        max = 1000,
        value = 500
      ),
      sliderInput(
        "width2",
        "Graph Width",
        min = 200,
        max = 1000,
        value = 500
      ),
      numericInput("xaxislowerlimits2", "X axis lower Limit", value = 600),
      numericInput("xaxislimits2", "X axis upper Limit", value = 1000),
      numericInput("xaxisinterval2", "X axis interval", value = 200),
      numericInput("yaxislimits2", "Y axis upper limit", value = 6000),
      numericInput("yaxisinterval2", "Y axis interval", value = 1000),
      
      div(style = "height:22px", h4("Filtered Out")),
      div(style = "display:inline-block;width:49%;height:60px;text-align:center", colourInput("mzSpeccolorchoices1", "Color", value = "white")),
      div(style = "height:22px", h4("No Overlap")),
      div(style = "display:inline-block;width:49%;height:60px;text-align:center", colourInput("mzSpeccolorchoices2", "Color", value = "green")),
      div(style = "height:22px", h4("1 Overlap")),
      div(style = "display:inline-block;width:49%;height:60px;text-align:center", colourInput("mzSpeccolorchoices3", "Color", value = "orange")),
      div(style = "height:22px", h4("2 Overlap")),
      div(style = "display:inline-block;width:49%;height:65px;text-align:center", colourInput("mzSpeccolorchoices4", "Color", value = "red")),
      div(style = "height:22px", h4("3 Overlap")),
      div(style = "display:inline-block;width:49%;height:65px;text-align:center", colourInput("mzSpeccolorchoices5", "Color", value = "purple")),
      
      div(style = "height:22px", h4("Background Color")),
      div(style = "display:inline-block;width:49%;text-align:center", colourInput("backgroundchoices2", "Background", value = "white"))
    ),
    mainPanel(
          actionButton("plotButton2", "Press to Plot or Apply Changes"),
      downloadButton("downloadbutton2", "Download Mass Spectra"),
      plotOutput(outputId = "plot2")
    ),
  ),
)))
 
#server####
server <- function(input, output, session) {
  
#server main         ####
  
  datafile <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    switch(
      ext,
      csv = vroom::vroom(input$upload$datapath, delim = ","),
      tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
    
  })
  
  a1 <- eventReactive(input$button1, {
    input$inputA
  })
  b1 <- eventReactive(input$button1, {
    input$inputB
  })
  c1 <- eventReactive(input$button1, {
    input$inputC
  })
  d1 <- eventReactive(input$button1, {
    input$intensityInputA
  })
  e1 <- eventReactive(input$button1, {
    input$intensityInputB
  })
  f1 <- eventReactive(input$button1, {
    input$RC_FilterF
  })
  g1 <- eventReactive(input$button1, {
    input$inputG
  })
  h1 <- eventReactive(input$button1, {
    input$inputH
  })
  
  a2 <- eventReactive(input$button2, {
    input$inputA
  })
  b2 <- eventReactive(input$button2, {
    input$inputB
  })
  c2 <- eventReactive(input$button2, {
    input$inputC
  })
  d2 <- eventReactive(input$button2, {
    input$intensityInputA
  })
  e2 <- eventReactive(input$button2, {
    input$intensityInputB
  })
  f2 <- eventReactive(input$button2, {
    input$RC_FilterF
  })
  g2 <- eventReactive(input$button1, {
    input$inputG
  })
  h2 <- eventReactive(input$button1, {
    input$inputH
  })
  #eventreactive function for buttons
  
  mz_RKMD_table <-
    reactive({
      arrangeKMD_Rcode(datafile(), a1(), b1(), c1(), d1(), e1(),f1(),g1(),h1())
    })
  RKMD_mz_table <-
    reactive({
      arrangeKMD_RKMDvsmz_count_repeats(datafile(), a2(), b2(), c2(), d2(), e2(),f2(),g2(),h2())
    })
  
  output$TableOutput <- renderDataTable({
    mz_RKMD_table()
  })
  output$TableOutput2 <- renderDataTable({
    RKMD_mz_table()
  })
  #make sure the function itself does NOT have read.csv since the file input already reads the csv file in
  
  KMD_Analysis_Template <- reactive({
    tibble(
      RKMD_Name = NA,
      RKMD_TheoMass = NA,
      blank1 = NA,
      blank2 = NA,
      blank3 = NA,
      blank4 = NA,
      blank5 = NA,
      MassList = NA,
      Intensity = NA
    )
  })
  
  output$download0 <- downloadHandler(
    filename = function() {
      paste("KMD_Analysis_Template", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(KMD_Analysis_Template(),
                na = "",
                row.names = FALSE,
                file)
    }
  )
  
  #downloading the KMD Analysis Template
  output$download1 <- downloadHandler(
    filename = function() {
      paste("calculated mz_vs_RKMD", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(mz_RKMD_table(), na = "", file)
    }
  )
  output$download2 <- downloadHandler(
    filename = function() {
      paste("calculated RKMD_vs_mz", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(RKMD_mz_table(), na = "", file)
    }
  )
#graphing server      ####
  plot1 <- eventReactive(input$plotButton, {
    generate(
      RKMD_mz_table(),
      shapelistchoices(),
      colorlistchoices(),
      sizechoice(),
      backgroundchoice(),
      linechoice(),
      inputxlower1(),
      inputxupper1(),
      inputylower1(),
      inputyupper1()
    )
  })
  
  plotlist <- reactive({
    plot1()[[1]]
  })
  titlelist <- reactive({
    plot1()[[2]]
  })
  plotchoices1 <- reactive({
    input$plotchoices
  })
  
  shapelistchoices <-
    reactive({
      as.integer(
        c(
          input$shapechoices1,
          input$shapechoices2,
          input$shapechoices3,
          input$shapechoices4
        )
      )
    })
  colorlistchoices <-
    reactive ({
      as.character(
        c(
          input$colorchoices1,
          input$colorchoices2,
          input$colorchoices3,
          input$colorchoices4
        )
      )
    })
  sizechoice <- reactive(input$sizechoices1)
  backgroundchoice <- reactive(input$backgroundchoices1)
  linechoice <- reactive(input$linechoices1)
  inputxlower1 <- reactive (input$inputxlower)
  inputxupper1 <- reactive (input$inputxupper)
  inputylower1 <- reactive (input$inputylower)
  inputyupper1 <- reactive (input$inputyupper)
  
  output$shapeKey <- renderPlot(show_point_shapes())
  output$plotselect <- renderUI({
    selectInput(inputId = "plotchoices", "Select RKMD", titlelist())
  })
  
  output$plot1ID <- renderPlot(
    width = function()
      input$width,
    height = function()
      input$height,
    res = 96,
    {
      plotlist()[[input$plotchoices]]
    }
  )
  output$downloadg2 <- downloadHandler(
    filename = function() {
      paste("RKMD_vs_mz_",
            input$plotchoices,
            "_",
            Sys.Date(),
            ".png",
            sep = "")
    },
    content = function(file) {
      png(
        file,
        width = input$width * 2,
        height = input$height * 2,
        res = 100 * input$sizechoices1
      )
      plot(plotlist()[[input$plotchoices]])
      dev.off()
    }
  )
  
#msSpec graphing server ####

colorlistchoices2 <-
  reactive ({
    as.character(
      c(
        "0"=input$mzSpeccolorchoices1,
        "1"=input$mzSpeccolorchoices2,
        "2"=input$mzSpeccolorchoices3,
        "3"=input$mzSpeccolorchoices4,
        "[4-10000000]"=input$mzSpeccolorchoices5
      )
    )
  })
backgroundchoice2 <- reactive(input$backgroundchoices2)
xaxislowerlimits2 <- reactive(input$xaxislowerlimits2)
xaxislimits2 <- reactive(input$xaxislimits2)
xaxisinterval2 <- reactive(input$xaxisinterval2)
yaxislimits2 <- reactive(input$yaxislimits2)
yaxisinterval2 <- reactive(input$yaxisinterval2)


plot2 <- eventReactive(input$plotButton2, {
  graphSpectrum(
    tibble(x=RKMD_mz_table()[[1]],y=datafile()[[9]],z=RKMD_mz_table()[[3]]),
    colorlistchoices2(),
    backgroundchoice2(),
    xaxislowerlimits2(),xaxislimits2(), xaxisinterval2(),
    yaxislimits2(), yaxisinterval2()
  )
})

output$plot2 <- renderPlot(
  width = function()
    input$width2,
  height = function()
    input$height2,
  res = 96,
  {
    plot2()
  }
)

output$downloadbutton2 <- downloadHandler(
  filename = function() {
    paste("mzSpec",
          "_",
          Sys.Date(),
          ".png",
          sep = "")
  },
  content = function(file) {
    png(
      file,
      width = input$width2 * 2,
      height = input$height2 * 2
    )
    plot(plot2())
    dev.off()
  }
)


}

shinyApp(ui = ui, server = server)











