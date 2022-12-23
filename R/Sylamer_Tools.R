#Ines Patop, 2022
#The things you can use for design of shRNAs

#function to get "not in" from %in% in dplyr
'%!in%' <- function(x,y)!('%in%'(x,y))

#' Sort differential gene expression resutls
#'
#' This function sort results from differential gene expression by any parameter.
#'
#'
#' @param DEres table with differential expression, one column must be gene names and must be named "gene", there has to be at least another column that will be used to sort genes
#' @param sortBy name of column used to sort, default is "log2FoldChange"
#' @param outName name of the file to write the sorted list to. Deafult is ="ForSylamer.txt"
#'
#' @param return Boolean, if True, then return the sorted list. Default is False.
#'
#' @return return the sorted list of genes list to run in sylamer
#'
#' @export
#'
#' @examples
#'
#' #sort_forSylamer(DEres="res.txt",sortBy="log2FoldChange",outName="ForSylamer.txt")
#'
#'
#'

sort_forSylamer <- function(DEres="res.txt",sortBy="log2FoldChange",outName="ForSylamer.txt",return=F){
  a<-read.delim(DEres)
  a<-a[,c("gene",sortBy)]
  a = a[order(a[,2], decreasing = F),]
  write.table(a, outName, row.names = F, sep = "\t", quote = F)
  if(return){return(a)}
}

#' Generate shell script to run Sylamer in one or multiple sorted lists of genes
#'
#' @param utrFile File with the UTR sequences. Defauls is "dm3_flybase_3utrs.fa"
#' @param kSylamer Oligo letter length or Word size, has to be the same length as those in the "words" file. Default is 6
#' @param kMarkov Size of a smaller word to be used for correcting composition biases. Default is 4
#' @param winSize How many sequences should be added in each consecutive window. Default is 200
#' @param sylamer Sylamer executable path. Default is "/opt/linux64bin/sylamer"
#' @param sylOutput Name of output file to write table in. Default is "shRNA_sylamer.output.tab"
#' @param imageOut Name of output file to write plot in. Default is "shRNA_image.output.pdf"
#' @param DEsorted Name of one or multiple list of sorted genes to input in sylamen. Default is c("ForSylamer.txt","ForSylamer_2.txt")
#' @param extras Any extra parametere to pass into Sylamer. Default is "-a 10 -aa 5 --funny-ok"
#' @param overwrite Boolean. If True, the script will overwrite any existing script with the same name. Default is T
#' @param sylamer_comand Name of the file to write the command into "sylamer_command.sh"
#' @return a file ready to run "sh sylamer_command.sh"
#'
#' @export
#'
#' @examples
#'
#' #create_sylamer_script(utrFile="dm3_flybase_3utrs.fa",kSylamer=6,kMarkov =4,winSize=200,sylamer="/opt/linux64bin/sylamer",sylOutput="shRNA_sylamer.output.tab",imageOut="shRNA_image.output.pdf",DEsorted=c("ForSylamer.txt","ForSylamer_2.txt"), extras = "-a 10 -aa 5 --funny-ok", overwrite=T,sylamer_comand="sylamer_command.sh")
#'
#'

create_sylamer_script <- function(utrFile="dm3_flybase_3utrs.fa",kSylamer=6,kMarkov =4,winSize=200,sylamer="/opt/linux64bin/sylamer",sylOutput="shRNA_sylamer.output.tab",imageOut="shRNA_image.output.pdf",DEsorted=c("ForSylamer.txt","ForSylamer_2.txt"), extras = "-a 10 -aa 5 --funny-ok", overwrite=T,sylamer_comand="sylamer_command.sh"){

  #setup error
  if(!file.exists(utrFile)){
    stop('UTR file not readable')
  }

  for(i in 1:length(DEsorted)){
    print(paste("Reading file: ",DEsorted[i]))
    if(!file.exists(DEsorted[i])){
      stop('DE file file not readable')
    }
  }

  #remove file
  if(overwrite){
    #rm(list = ls(pattern = "ForSylamer"))
    file.remove(sylamer_comand)
  }

  for(i in DEsorted){
      cat(paste0(sylamer," -fasta ",utrFile," -grow ",winSize," -k ",kSylamer," --print-id -m ", kMarkov," -subset ",i," -o ",sylOutput),file = sylamer_comand,append = T,sep = "\t",fill = T)
    }

}

#' Plot the output of Sylamer
#'
#' @param oligosChosen Sequences to highlight (Any particular words that you want highlighted) Defaults is =c("TGTAAA","GATGCT")
#' @param sylamer.out.tab Name of sylamer output to read from. Dafault is "shRNA_sylamer.output.tab"
#' @param topOligos Number of top oligos to plot. How many best words should I highlight in the image? Defauls is 0
#' @param lowOligos Number of low oligos to plot. Defauls is 0
#' @param outName name for the plot output. Default is "sylamerplot"
#' @param greylines Boolean, To plot or not greay background oligos. Default is F
#'@param kSylamer Oligo letter length or Word size, has to be the same length as those in the "words" file. Default is 6
#' @param kMarkov Size of a smaller word to be used for correcting composition biases. Default is 4
#' @param free_ylim Boolean, if TRUE then y lim will be the min and max values from Sylamr output. If FALSE, then ylim=(-10,10)
#'
#' @return a file ready to run "sh sylamer_command.sh"
#'
#' @export
#'
#' @examples
#'
#' #create_sylamer_script(utrFile="dm3_flybase_3utrs.fa",kSylamer=6,kMarkov =4,winSize=200,sylamer="/opt/linux64bin/sylamer",sylOutput="shRNA_sylamer.output.tab",imageOut="shRNA_image.output.pdf",DEsorted=c("ForSylamer.txt","ForSylamer_2.txt"), extras = "-a 10 -aa 5 --funny-ok", overwrite=T,sylamer_comand="sylamer_command.sh")
#'
#'

sylamerplots<-function(oligosChosen=c("TGTAAA","GATGCT"),sylamer.out.tab="shRNA_sylamer.output.tab",topOligos=0,lowOligos=0,outName="sylamerplot",greylines=F,kSylamer=6,kMarkov =4,winSize=200,plot=T,free_ylim=F){

  # How many best words should I highlight in the image?
  topOligos <- topOligos
  lowOligos <- lowOligos

  # Any particular words that you want highlighted
  chosenOligos<- oligosChosen

  #read reults table
  sylTable<- read.table(file = sylamer.out.tab , sep="\t", row.names=1, header=T, check.names=F)

  xVals<- as.numeric(colnames(sylTable))

  sylTable <- cbind("0"=0,sylTable)   # To add an initial column of 0s

  xVals<- as.numeric(colnames(sylTable))

  # Generate the plot area
  yMin   <- min(sylTable)
  yMax   <- max(sylTable)
  yRange <- yMax - yMin

  pdf(paste0(outName,".pdf"))


  if(free_ylim){
    .ylim=c(yMin,yMax)
  }else{
    .ylim=c(-10,10)
  }
  plot(NULL, xlab="3UTR sequences sorted", ylab="log10(enrichment P-value)", axes=T,  main=paste0("Sylamer landscape using words of length: ",kSylamer, " for ",outName),ylim=.ylim,xlim=range(xVals) )
  #ylim=c(round(yMin-yRange/10),round(yMax+yRange/10)), xlim=range(xVals))

  # It can save time to plot no more than ~1,000 lines (particularly for all words of length 7 or 8)
  if (greylines == T){
    maxPlot <- nrow(sylTable)
    # Plot the background lines
    for (i in 1:maxPlot) {
      lines(xVals,sylTable[i,], col='grey')
    }

  } else{
    maxPlot = 0
  }


  # Draw a reference line at 0
  abline(h=0)
  abline(h=2.3, lty=2)
  abline(h=-2.3, lty=2)

  # Up/Down best words
  oligosUp <- c()
  oligosDown <- c()
  oligosChosen <- c()
  if (topOligos > 0) { # Only if I really want these plots
    oligosUp <- names((sort(apply(sylTable,1, function(x) {max(x[is.finite(x)])}),decreasing=TRUE))[1:topOligos])
  }
  if (lowOligos > 0) {
    oligosDown <- names((sort(apply(sylTable,1,function(x) {min(x[is.finite(x)])}),decreasing=FALSE))[1:lowOligos])
  }
  if (length(chosenOligos) > 0) {
    oligosChosen <- chosenOligos[chosenOligos %in% rownames(sylTable)]
  }
  oligosAll <- unique(c(oligosUp,oligosDown,oligosChosen))
  oligosAll <- names(sort(apply(sylTable,1, function(x) {max(abs(x[is.finite(x)]))})[oligosAll],decreasing=TRUE))

  if (length(oligosAll) > 0) {
    colors   <- rainbow(length(oligosAll))
    names(colors) <- oligosAll
    for (i in rev(seq_along(oligosDown))) {
      lines(xVals, sylTable[oligosDown[i],], col=colors[oligosDown[i]], lwd=2)
    }
    for (i in rev(seq_along(oligosUp))) {
      lines(xVals, sylTable[oligosUp[i],], col=colors[oligosUp[i]], lwd=2)
    }
    for (i in rev(seq_along(oligosChosen))) {
      lines(xVals, sylTable[oligosChosen[i],], col=colors[oligosChosen[i]], lwd=2)
    }
    if (topOligos >0) {
      legend('topleft', inset=c(0.01,0.01), legend=oligosUp, lwd=2, lty=1, horiz=TRUE, col=colors[oligosUp],
             cex=0.6, bg='white', title="Words with highest peak")
    }
    if (lowOligos >0) {
      legend('bottomleft', inset=c(0.01,0.01), legend=oligosDown, lwd=2, lty=1, horiz=TRUE, col=colors[oligosDown],
             cex=0.6, bg='white', title="Words with lowest peak")
    }
    if (length(oligosChosen) >0) {
      legend('topright', inset=c(0.01,0.01), legend=oligosChosen, lwd=2, lty=1, horiz=F, ncol=1, col=colors[oligosChosen],
             cex=0.6, bg='white', title="Selected words")
    }
  }
  dev.off()

  if(plot==T){

    if(free_ylim){
      .ylim=c(yMin,yMax)
    }else{
      .ylim=c(-10,10)
    }

    plot(NULL, xlab="3UTR sequences sorted", ylab="log10(enrichment P-value)", axes=T,  main=paste0("Sylamer landscape using words of length: ",kSylamer, " for ",outName),,ylim=.ylim,xlim=range(xVals) )


    #ylim=c(round(yMin-yRange/10),round(yMax+yRange/10)), xlim=range(xVals))

    # It can save time to plot no more than ~1,000 lines (particularly for all words of length 7 or 8)
    if (greylines == T){
      maxPlot <- nrow(sylTable)
      # Plot the background lines
      for (i in 1:maxPlot) {
        lines(xVals,sylTable[i,], col='grey')
      }

    } else{
      maxPlot = 0
    }


    # Draw a reference line at 0
    abline(h=0)
    abline(h=2.3, lty=2)
    abline(h=-2.3, lty=2)

    # Up/Down best words
    oligosUp <- c()
    oligosDown <- c()
    oligosChosen <- c()
    if (topOligos > 0) { # Only if I really want these plots
      oligosUp <- names((sort(apply(sylTable,1, function(x) {max(x[is.finite(x)])}),decreasing=TRUE))[1:topOligos])
    }
    if (lowOligos > 0) {
      oligosDown <- names((sort(apply(sylTable,1,function(x) {min(x[is.finite(x)])}),decreasing=FALSE))[1:lowOligos])
    }
    if (length(chosenOligos) > 0) {
      oligosChosen <- chosenOligos[chosenOligos %in% rownames(sylTable)]
    }
    oligosAll <- unique(c(oligosUp,oligosDown,oligosChosen))
    oligosAll <- names(sort(apply(sylTable,1, function(x) {max(abs(x[is.finite(x)]))})[oligosAll],decreasing=TRUE))

    if (length(oligosAll) > 0) {
      colors   <- rainbow(length(oligosAll))
      names(colors) <- oligosAll
      for (i in rev(seq_along(oligosDown))) {
        lines(xVals, sylTable[oligosDown[i],], col=colors[oligosDown[i]], lwd=2)
      }
      for (i in rev(seq_along(oligosUp))) {
        lines(xVals, sylTable[oligosUp[i],], col=colors[oligosUp[i]], lwd=2)
      }
      for (i in rev(seq_along(oligosChosen))) {
        lines(xVals, sylTable[oligosChosen[i],], col=colors[oligosChosen[i]], lwd=2)
      }
      if (topOligos >0) {
        legend('topleft', inset=c(0.01,0.01), legend=oligosUp, lwd=2, lty=1, horiz=TRUE, col=colors[oligosUp],
               cex=0.6, bg='white', title="Words with highest peak")
      }
      if (lowOligos >0) {
        legend('bottomleft', inset=c(0.01,0.01), legend=oligosDown, lwd=2, lty=1, horiz=TRUE, col=colors[oligosDown],
               cex=0.6, bg='white', title="Words with lowest peak")
      }
      if (length(oligosChosen) >0) {
        legend('topright', inset=c(0.01,0.01), legend=oligosChosen, lwd=2, lty=1, horiz=F, ncol=1, col=colors[oligosChosen],
               cex=0.6, bg='white', title="Selected words")
      }
    }
  }
}


