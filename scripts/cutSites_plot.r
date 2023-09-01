library(ggplot2)
library(dplyr)
library(plyr)
args <- commandArgs(trailingOnly = TRUE)
ID = args[1]


# combine for ALU
Forward = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGA"
Forward.list <- list()
Reverse.list <- list()
for(start in 0:15){
    file = sprintf("read_start%s_cutsite_Alu.txt", start)
    file2 = sprintf("read_start%s_cutsite_Alu_for_plot.txt", start)
    system(paste0("grep Alu ", file, " > ", file2))
    site.tb <- read.table(file2, blank.lines.skip=T,sep = '\t', fill=T, stringsAsFactors = F)
    site.tb <- site.tb[site.tb$V4=="Alu",]
    Forward.tb <- site.tb[grepl("Forward",site.tb$V2),]
    Reverse.tb <- site.tb[grepl("Reverse",site.tb$V2),]
    Forward.list[[as.character(start)]] <- Forward.tb[, c("V1", "V3")]
    Reverse.list[[as.character(start)]] <- Reverse.tb[, c("V1", "V3")]
}  
# length(Reduce(union, lapply(Forward.list, function(x){x$V1})))
# length(Reduce(union, lapply(Reverse.list, function(x){x$V1})))
test <- join_all(Forward.list, type="full")
# tail(sort(table(test2$V1)))
test <- test[order(test$V1, test$V3), ]
test <- test[!duplicated(test$V1), ]

test2 <- join_all(Reverse.list, type="full")
# tail(sort(table(test2$V1)))
test2 <- test2[order(test2$V1, test2$V3, decreasing=T), ]
test2 <- test2[!duplicated(test2$V1), ]

# Forward
Forward.pp <- data.frame(site=1:60, number=0)
t <- table(test$V3)[names(table(test$V3)) %in% 1:60]
Forward.pp[names(t),"number"] <- t
Forward.pp$site <- factor(Forward.pp$site)
Forward.pp$number <- as.numeric(Forward.pp$number)
# reverse
Reverse.pp <- data.frame(site=1:60, number=0)
t <- table(test2$V3)[names(table(test2$V3)) %in% 1:60]
Reverse.pp[names(t),"number"] <- t
Reverse.pp$site <- factor(Reverse.pp$site)
Reverse.pp$number <- as.numeric(Reverse.pp$number)


#
pp <- rbind(data.frame(Forward.pp, Orientation="Forward"),
          data.frame(Reverse.pp, Orientation="Reverse"))
pp[pp$Orientation=="Reverse", "number"] <- -pp[pp$Orientation=="Reverse", "number"]


ggplot(pp, aes(x=site, y=number, fill=Orientation)) +
    geom_bar(stat="identity") +
    # facet_wrap(~class, nrow = 2, scales="free") +
    scale_fill_manual(values = c("#4678BA", "#EF8354")) +
    ggtitle(ID) +
    theme_bw() +
    theme(
      panel.grid=element_blank(),
      axis.text.x = element_text(angle = 0, color = "black",
                                 size = 8, hjust = 0.5, vjust=0.5),
      plot.title=element_text(size=14, face="bold")
    ) +
    scale_x_discrete(labels=strsplit(Forward, "")[[1]])+
    ylab("the number of reads")


ggsave(sprintf("%s.Alu_cleavage_sites.pdf", ID), width=8, height=5)



