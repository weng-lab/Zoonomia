
# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making wordCloud from GO analysis results.

library(ggplot2)

library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")

args = commandArgs(trailingOnly=TRUE)
GO_result = args[1]
prefix = args[2]

####################
# ubi-PLS vs. PLS
####################
dat = read.table(GO_result,sep="\t", header = TRUE)

#--------- plotting enriched terms
d1 = dat[dat[,5]=="+",]
d1_q = d1[d1[,8]<1e-2,] # adjusted enriched p-value cutoff

## get terms
text = as.vector(d1_q[,1])

## clean the text
docs <- Corpus(VectorSource(text))
inspect(docs)
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, ":")
docs <- tm_map(docs, removeWords, c("process"))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, removeNumbers)

## Build a term-document matrix
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)

# set.seed(1234)
pdf(paste(prefix,"_pos.pdf",sep=""))
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=100, random.order=FALSE, rot.per=0.35,
          colors="black")
dev.off()

#--------- plotting depleted terms
d2 = dat[dat[,5]=="-",]
d2_q = d2[d2[,8]<1e-5,]

## get terms
text = as.vector(d2_q[,1])

## clean the text
docs <- Corpus(VectorSource(text))
inspect(docs)
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, ":")
docs <- tm_map(docs, removeWords, c("process"))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, removeNumbers)


## Build a term-document matrix
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)

# set.seed(1234)
pdf(paste(prefix,"_neg.pdf",sep=""))
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=100, random.order=FALSE, rot.per=0.35,
          colors="black")
dev.off()
