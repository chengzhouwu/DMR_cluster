# the programming is used to cluster the CpG sites

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(geograbi) 
library(limma)
library(corrplot)

# functions used in the programming 
#input: x, y are both n by p matrices 
#output: a vector of length p indicating the column-wise correlation between x and y.

colCors = function(x, y) {
        sqr = function(x) x*x
        if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
                stop("Please supply two matrices of equal size.")
        x   = sweep(x, 2, colMeans(x))
        y   = sweep(y, 2, colMeans(y))
        cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
        return(cor)
}




load('beta.Rdata')
t.beta <- t(beta) 

#short the long row name of t.beta to identify with samples and vars
rownames(t.beta) <- substring(rownames(t.beta),1,10)
colnames(beta) <- substring(colnames(beta),1,10)


#get the map information for CpG sites 
map.data.subset <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
colnames(map.data.subset)
#delete X and Y since our data delete the chr X and Chr Y 
map.data.subset = map.data.subset[map.data.subset$chr!='chrY'&map.data.subset$chr!='chrX',] 
map.data.subset <-map.data.subset$chr[map.data.subset$Name %in% rownames(beta)]
table(map.data.subset)

Z1 <- t.beta
dim(Z1)

################################
# get the block information ###
#############################
##get start.loc.B and end.loc.B
##get start.loc.B and end.loc.B
start.loc.B=c()
end.loc.B=c()
count=0
for (chr.num in names(table((map.data.subset)))){
        count=count+1
        tbl=table((map.data.subset))
        
        which(names(tbl)==chr.num)
        index=(cumsum(tbl)[count]-tbl[count]+1):cumsum(tbl)[count]
        
        
        record_cor=colCors(Z1[,index[-length(index)]],Z1[,index[-1]])
        
        
        consecutive.table=rle(as.integer(record_cor>0.5))
        block.len=consecutive.table[[1]]
        block.val=consecutive.table[[2]]
        
        temp.loc=1:length(index)
        temp.block.start=rep(0,sum(block.val==1,na.rm=T))
        temp.block.end=rep(0,sum(block.val==1,na.rm=T))
        for (i in 1:sum(block.val==1,na.rm = T)){
                temp.block.start[i]=(cumsum(block.len)[block.val==1]-block.len[block.val==1]+1)[i]
                temp.block.end[i]=cumsum(block.len)[block.val==1][i]+1
                
        }
        
        temp.start.loc.B=1:length(index)
        temp.end.loc.B=1:length(index)
        if (sum(block.val==1,na.rm = T)>0){
                for (i in 1:sum(block.val==1,na.rm = T)){
                        temp.start.loc.B[temp.block.start[i]:temp.block.end[i]]=temp.block.start[i]
                        temp.end.loc.B[temp.block.start[i]:temp.block.end[i]]=temp.block.end[i]
                }
                temp.start.loc.B=unique(temp.start.loc.B)
                temp.end.loc.B=unique(temp.end.loc.B)
        }
        
        
        if (count==1){
                start.loc.B=temp.start.loc.B
                end.loc.B=temp.end.loc.B
        }else{
                start.loc.B=c(start.loc.B,temp.start.loc.B+cumsum(tbl)[count-1])
                end.loc.B=c(end.loc.B,temp.end.loc.B+cumsum(tbl)[count-1])
        }
        
}

table(end.loc.B-start.loc.B+1)
total_compare_LR <- sum(table(end.loc.B-start.loc.B+1))

L.vec <- as.numeric(names(table(end.loc.B-start.loc.B+1)))
length(L.vec)



