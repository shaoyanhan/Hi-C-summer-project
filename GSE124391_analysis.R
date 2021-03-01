#初始化文件信息
path<-"/store/yhshao/GSE124391_RAW";#交互信息文件的路径
resolution<-100000#全局分辨率
fileNames<-dir(path) ;#按照路径获取文件名
fileNames<-fileNames[-c(1:5)]#去除前五个（多细胞Hi-C数据以及其他不相关文件）
sana_mESC<-read.csv(file = '/store/yhshao/mESC.txt',header = T)#读入样本详细信息文件
sana_mESC<-sana_mESC[order(sana_mESC$GEO_Accession..exp.),]#对样本信息文件按照GEO号进行排序，因为路径里面的样本是按照GEO号排序的
sana_mESC<-sana_mESC[-c(1:3),]#去除样本信息文件中多细胞Hi-C样本的信息
fileNames<-fileNames[-which(sana_mESC$Organism=='Homo sapiens')]#将样本列表和样本信息文件中的人类样本去除
sana_mESC<-sana_mESC[-which(sana_mESC$Organism=='Homo sapiens'),]


#提取细胞类型名称
qc<-array()#qc是一个数组，用于存储每个样本的染色体内交互占比，正常情况下一个单细胞数据染色体内交互占比应该为50～60%
# 以下注释部分为计算占比方法，由于较为耗时，只建议跑一次然后使用结果文件（或者直接使用下面我路径中已跑好的文件）
# j=1
# for (filenum in 1:length(filePath)) {#遍历所有样本文件
#   data1<-0#初始化存储样本文件的表格
#   data1<-read.table(filePath[filenum],header = F)#读入样本文件
#   colnames(data1)<-c('chr1','pos1','chr2','pos2')#写入各列名称
#   data1$chr1<-as.character(data1$chr1)#将第一号和第二号染色体列的数据类型转换为字符串类型
#   data1$chr2<-as.character(data1$chr2)
#   qc[j]<-length(which(data1$chr1==data1$chr2))/nrow(data1)#计算染色体内交互占比
#   j=j+1
#   print(paste0('第',filenum,'条样本统计完毕'))
# }
qc<-read.csv(file='/store/yhshao/Rwork/qc.csv')#读取占比文件
qc<-qc[,-1]#去除序号列
fileNames<-fileNames[-which(qc<0.5)]#将染色体内交互占比小于50%的样本去除
sana_mESC<-sana_mESC[-which(qc<0.5),]
groups_qc<-sana_mESC$source_name#将样本信息中每个样本的名称（即所属类别）提取出来
filePath <- sapply(fileNames, function(x){ paste(path,x,sep='/')})#生成样本路径列表
print(paste("共",length(filePath),"条待处理样本"))


#获取文件列表里面每一个样本文件的大小信息，以用于之后的矩阵归一化，消除文件大小对降维的影响
fs<-matrix(nrow = length(filePath),ncol = 2)#初始化文件大小信息矩阵：358*2
fs[,1]<-round((file.info(filePath)$size)/1024^2)#文件大小化至兆级并四舍五入
fs[,2]<-as.character(sana_mESC$source_name)#写入对应样本名称


#注：43-147行为compartment、tad、matsize_com、matsize_tad代码的复用
#读入分compartment结果文件并预处理
compartments<-0
compartments<-read.csv(file = '/store/yhshao/work/compartment/ES_E14_GRCm38_100kb_com.cis.vecs.tsv',header = T)
aa<-array()#aa用于记录文件中的空行
j=1
for (i in 1:nrow(compartments)) {#先将文件的每一行按照制表符拆开并判断是否为空行
  if(strsplit(as.character(compartments[i,1]),split = '\t')[[1]][4]==''){
    aa[j]<-i;j=j+1;
  }
}
compartments<-compartments[-aa,1]#将空行去除
bb<-matrix(nrow = length(compartments),ncol = 4)#生成一个4列的表格
for (i in 1:nrow(bb)) {
  bb[i,]<-strsplit(as.character(compartments[i]),split = '\t')[[1]][1:4]#将基准文件里每行数据写入表格
}
bb<-as.data.frame(bb)#将表格转化为数据框方便后续取值操作
bb$V2<-round(as.numeric(as.character(bb$V2))/resolution)#将compartment起始位点除以分辨率并四舍五入
bb$V2[which(bb$V2==0)]<-1#若起始位点为零，需要改为1方便之后计算
bb$V3<-round(as.numeric(as.character(bb$V3))/resolution)#将compartment结束位点除以分辨率并四舍五入
bb$V4<-as.numeric(as.character(bb$V4))#将compartment值转换为数值类型
colnames(bb)<-c('chrom','start','end','Evalue')#经过上述处理得到了一个表格，存储了染色体名称、compartment起始与结束位点、compartment值
compartments<-bb#表格改名
compartments$Evalue[which(compartments$Evalue>0)]<-1#将compartment值转换为正负1表示
compartments$Evalue[which(compartments$Evalue<0)]<--1
#读入分tad结果文件并预处理
tads<-0
tads<-read.csv(file = "/store/yhshao/work/tad/ES_E14_GRCm38_100kb_tad.txt",header = T)
aa<-array()
j=1#进行去除nan行操作，即将tad的区间合并，将所有有数值的行保留，则两两行之间即代表一个tad，下面的字符串拆分和表格填入操作与compartment类似
for (i in 1:nrow(tads)) {
  if(strsplit(as.character(tads[i,1]),split = '\t')[[1]][7]=='nan'){
    aa[j]<-i;j=j+1;
  }
}
tads<-tads[-aa,1]
bb<-matrix(nrow = length(tads),ncol = 3)
for (i in 1:nrow(bb)) {
  bb[i,]<-strsplit(as.character(tads[i]),split = '\t')[[1]][1:3]
}
bb<-as.data.frame(bb)
bb$V2<-round(as.numeric(as.character(bb$V2))/resolution)
bb$V2[which(bb$V2==0)]<-1
bb$V3<-round(as.numeric(as.character(bb$V3))/resolution)
colnames(bb)<-c('chrom','start','end')#同上，这里生成了一个表格用于存储染色体名称、起始与结束位点
tads<-bb

#计算热图矩阵的尺度，之后填充的热图矩阵需要基于基准文件的尺度，因此需要单独计算以达到规范性，下面的matsize_com与matsize_tad是21*4的表格
chrs<-seq(1,19)
chrs[20:21]<-c('X','Y')
chrs<-paste('chr',chrs,sep ='')#生成染色体名称序列
matsize_com<-matrix(nrow = length(chrs),ncol = 4)#compartment的热图矩阵尺度表格，这里一共4列，分别记录染色体名称、热图矩阵最大和最小值以及该条染色体上的compartment数目
i=1
a1=0
a2=0
for (chrom in chrs) {#向matsize_com填入数据
  if(length(which(compartments$chrom==chrom))==0){#如果基准文件里没有该染色体，将该染色体的热图矩阵尺度大小以及compartment数目记为-1
    matsize_com[i,1]<-chrom
    matsize_com[i,2:4]<-seq(-1,-1,length.out = 3)
    i=i+1
  }
  else{
    a1<-which(compartments$chrom==chrom)[1]#a1和a2记录了该条染色体在compartment基准文件里面的起始行和终止行
    a2<-which(compartments$chrom==chrom)[length(which(compartments$chrom==chrom))]
    s1<-compartments$start[a1]#s1和s2记录了该条染色体在compartment基准文件里面的起始位点和终止位点
    s2<-compartments$end[a2]
    j=0#用于记录compartment数目
    for (m in a1:a2) {#对该条染色体所拥有的compartment进行计数
      if(m==a2){#如果遍历到最后一行，则只有两种情况，要么最后一个compartment恰好只有一行，要么这是最后一个compartment的最后一行，所以直接加1即可
        j=j+1
        break()
      }
      if(compartments$Evalue[m]*compartments$Evalue[m+1]<0){#如果相邻两行的compartment值（+1/-1）异号，则意味着compartment的切换，也是加1
        j=j+1
      }
    }
    matsize_com[i,1]<-chrom#将染色体名称、热图矩阵尺度和compartment数目写入matsize_com
    matsize_com[i,2]<-s1
    matsize_com[i,3]<-s2
    matsize_com[i,4]<-j
    i=i+1
  }
}

matsize_tad<-matrix(nrow = length(chrs),ncol = 4)#tad的热图矩阵尺度表格，这里一共4列，分别记录染色体名称、热图矩阵最大和最小值以及该条染色体上的tad数目
i=1
a1=0
a2=0
for (chrom in chrs) {#向matsize_tad填入数据
  if(length(which(tads$chrom==chrom))<2){#因为tad基准文件为两行记录一个tad区间的格式，因此如果tad基准文件小于2，则表明基准文件在该分辨率下没有tad或者没有该染色体，在该行写入-1
    matsize_tad[i,1]<-chrom
    matsize_tad[i,2:4]<-seq(-1,-1,length.out = 3)
    i=i+1
  }
  else{
    a1<-which(tads$chrom==chrom)[1]#a1和a2记录了该条染色体在tad基准文件里面的起始行和终止行
    a2<-which(tads$chrom==chrom)[length(which(tads$chrom==chrom))]
    s1<-tads$start[a1]#s1和s2记录了该条染色体在tad基准文件里面的起始位点和终止位点
    s2<-tads$end[a2]
    matsize_tad[i,1]<-chrom#将染色体名称、热图矩阵尺度和tad数目写入matsize_tad
    matsize_tad[i,2]<-s1
    matsize_tad[i,3]<-s2
    matsize_tad[i,4]<-length(which(tads$chrom==chrom))-1#因为tad基准文件为两行记录一个tad区间的格式，因此每条染色体所具有的tad即为行数减1
    i=i+1
  }
}


#将计数矩阵每一行除以该样本文件的大小，达到归一化，这里使用tad矩阵示例
samplematrix_tad<-samplematrix_100k_qc_ES_mESC_tad
for (i in 1:nrow(fs)) {
  samplematrix_tad[i,]<-round(samplematrix_tad[i,]/as.numeric(fs[i,1]))#每行除以对应样本文件大小
}
zerocol<-0#该变量用于存储矩阵中数值全部为零的列号
j=1
for (i in 1:ncol(samplematrix_tad)) {
  if(sum(samplematrix_tad[,i])==0){
    zerocol[j]<-i
    j=j+1
  }
}
zerocol
#如果存在全部为零的列，要将其去除，因为零列无法进行降维
if(sum(zerocol)>0){
  samplematrix_tad<-samplematrix_tad[,-zerocol]
}


#下载可视化所需的包
install.packages("devtools",repo="http://cran.us.r-project.org")
install.packages("ggplot2")
library(devtools)
library(ggplot2)
install_github("vqv/ggbiplot")
library(ggbiplot)
install.packages("mnormt")
library(mnormt)
install.packages("psych")
library(psych)
install.packages("zip")
library(zip)
install.packages("factoextra")
library(factoextra)
install.packages("corrplot")
library("corrplot")
install.packages("Rtsne")
library(Rtsne)


#k-means聚类（效果较差，不推荐使用）
sample.pca<-principal(samplematrix,scores = T,nfactors = 2)
fviz_nbclust(sample.pca$scores,kmeans,method="wss")+geom_vline (xintercept=4,linetype=2)
pca_kmeans<-kmeans(sample.pca$scores,4)
fviz_cluster(pca_kmeans,data = sample.pca$scores)


#PCA可视化
sample.pca<-prcomp(samplematrix,scale. = T)#首先利用PCA计算函数之一“prcomp”对计数矩阵降维，结果返回一个list

fviz_eig(sample.pca,addlabels = TRUE) #碎石图,展示在不同PC上的方差解释度

# fviz_pca是factoextra包中的PCA可视化利器，其对“贡献度”的分析结果较为重要
# coord是坐标，与cor（相关性）数值相同，coord=loading * stdev=loadings*sqrt(eig)
# cos2为cor的平方，代表不同主成分对变量的代表性强弱，对特定变量，所有主成份上的cos2之和为1，其结果与cor类似
# 对向的是负相关。箭头越远离远原点、越靠近圆圈表明PC对其的代表性高（相关性强）

# 变量层面（ind）的代表性（cos2）和贡献度（contrib）
# axes:选择维度（PC）；col.var:选择cos2（变量）或ind（样本）；gradient.cols：颜色梯度；repel：利用算法使图中的文字不重合
# 按照cos2大小设定颜色梯度，也可以设置alpha梯度，在下面这个相关图中，靠近的变量表示正相关，直角为正交，大于直角为负相关；
fviz_pca_var(sample.pca,axes=c(1,2),col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE) 

# cos2的水滴图，行为变量，列为PC
corrplot(get_pca_var(sample.pca)$cos2, is.corr=FALSE)
# cos2在主成分上的加和，并排序，显示排名前10的结果
fviz_cos2(sample.pca, choice = "var", axes = 1:2 ,top = 10) 
# 下式可用于寻找最大最小结果，当变量较多时很方便
# which.max(get_pca_var(sample.pca)$contrib[,1]+get_pca_var(sample.pca)$contrib[,2])
# which.min(get_pca_var(sample.pca)$contrib[,1]+get_pca_var(sample.pca)$contrib[,2])

# contrib可视化：
# contrib是每个变量对某一特定PC的贡献，contrib=(var.cos2 * 100) / (total cos2 of the PC component)，
# 多个变量的contrib = [(Contrb1 * Eig1) + (Contrib2 * Eig2)]/(Eig1 + Eig2)
corrplot(get_pca_var(sample.pca)$contrib, is.corr=FALSE)
fviz_contrib(sample.pca, choice = "var", axes = 1:2 ,top = 10)

# 样本层面（ind）的代表性（cos2）和贡献度（contrib）
fviz_cos2(sample.pca, choice = "ind",axes=1:2)
fviz_contrib(sample.pca, choice = "ind", axes = 1:2)

# 普通的pca样本点聚类图
sample.pca<-prcomp(samplematrix,scale. = T)
# 利用ggbiplot进行可视化，groups是计数矩阵里面每一行样本的名称（即类别），choice：选取的PC
ggbiplot(sample.pca, groups = groups_qc,obs.scale = 1, var.scale = 1,var.axes = F,choices = c(1,2))

# 较为详尽的pca样本点聚类图
# 利用fviz_pca的结果，将样本点的大小代表cos2值，并将样本名称、类别显示出来,geom:添加样本名称（fill.ind），
# addEllipses：添加聚类的椭圆，col.ind：线颜色，pointsize：点大小，更多参数可以查看帮助?fviz_pca_ind
fviz_pca_ind(sample.pca, geom=c("point","text"),addEllipses = T,pointshape=21,col.ind="black",
             pointsize="cos2",fill.ind = groups_qc,repel = TRUE,axes = c(2,3))

# t-SNE可视化
# dims维度大小，max_iter迭代次数，perplexity参数的取值必须小于(nrow(data) - 1 )/ 3; 
# theta参数取值越大，结果的准确度越低，默认值为0.5；pca参数表示是否对输入的原始数据进行PCA分析，
# 然后使用PCA得到的topN主成分进行后续分析，t-SNE算法的计算量是特别大的，对于维度较高的数据数据，
# 先采用PCA降维可以有效提高运行的效率，默认采用top50的主成分进行后续分析，当然也可以通过initial_dims参数修改这个值
tsne_out <-Rtsne(samplematrix,dims=2,max_iter=2000,perplexity=40,pca=FALSE, theta=0.0)
ggplot(data=as.data.frame(tsne_out$Y),mapping=aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2],color=(groups_qc),shape=(groups_qc)))+geom_point()
#library(ggrepel);
#ggplot(data=as.data.frame(tsne_out$Y),mapping=aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2],color=(groups_qc),shape=(groups_qc)))+geom_point()+geom_text_repel(aes( label = rownames(as.data.frame(tsne_out$Y))),colour="black",size=2,alpha=0.5,segment.size=0.1,segment.alpha=0.5)


#tad或者染色体层面的热图绘制
filenum<-46#指定绘图的样本序号
filesize<-as.numeric(fs[filenum,1])#记录该样本文件的大小
chrom<-'chr8'#指定需要绘制热图的染色体
startsite<-362#指定绘制热图的区间
endsite<-372
data1<-0#初始化文件存储变量
data1<-read.table(filePath[filenum],header = F)#读入交互信息文件
colnames(data1)<-c('chr1','pos1','chr2','pos2')#写入列名称
data1$chr1<-as.character(data1$chr1)#将染色体名称列转化为字符串类型
data1$chr2<-as.character(data1$chr2)
data1$interval1<-round(data1$pos1/resolution)#将交互位点数值除以分辨率并四舍五入，存入interval列
data1$interval2<-round(data1$pos2/resolution)
s1=0
s2=0
a1=0
a2=0
s1<-as.numeric(matsize_tad[which(matsize_tad[,1]==chrom),2])#s1和s2用于存储绘制该染色体热图的矩阵大小（注：热图的绘制尺度是基于基准文件的）
s2<-as.numeric(matsize_tad[which(matsize_tad[,1]==chrom),3])
y<-matrix(0,nrow =s2,ncol =s2)#初始化热图矩阵
a1<-which(data1$chr1==chrom)#找出该样本文件中第一号染色体名称对应我们选取的染色体名称的行
a2<-which(data1$chr2[a1]!=chrom)#在上述行（a1）中找出第二号染色体名称不等于第一号染色体名称的行
if(length(a2)>0){
  data1<-data1[-a1[a2],]#若a2不为零，那么代表存在染色体间交互数据，将其去除
}
if(length(which(data1$interval1[which(data1$chr1==chrom)]>s2))>0){#找到样本文件中第一个交互位点的尺度超过基准文件尺度的行并去除
  data1<-data1[-which(data1$chr1==chrom)[which(data1$interval1[which(data1$chr1==chrom)]>s2)],]
}
if(length(which(data1$interval2[which(data1$chr1==chrom)]>s2))>0){#找到样本文件中第二个交互位点的尺度超过基准文件尺度的行并去除
  data1<-data1[-which(data1$chr1==chrom)[which(data1$interval2[which(data1$chr1==chrom)]>s2)],]
}
for (i in which(data1$chr1==chrom)) #使用文件中剩余的的染色体内交互填充热图矩阵
{
  y[data1[i,5],data1[i,6]]<-1+y[data1[i,5],data1[i,6]];
}
for (i in which(data1$chr1==chrom)) 
{
  y[data1[i,6],data1[i,5]]<-1+y[data1[i,6],data1[i,5]];
}
y<-log((y/filesize)+1)#将热图矩阵除以该样本文件的大小进行归一化，然后取log（拉开差异）并加1（为零的点无法取log）
#y[which(y>0)[which.min(y[y>0])]]#这两行代码可以查看矩阵中的非零最小值与最大值
#y[which.max(y)]

#绘制热图
ramp <- colorRamp(c("red", "white"));
cols<-rgb( ramp(seq(0, 1, length = 128)), max = 255);
tit<-paste('sample',filenum,'_',chrom,'(',startsite,'-',endsite,'/',resolution/1000,'k)',sep = '')#编辑标题
heatmap(y[startsite:endsite,startsite:endsite],col=rev(cols), Rowv=NA,Colv=NA, labRow=NA,labCol=NA,scale="none",main=tit);


#compartment层面的热图
#下面这段代码是填充热图矩阵，与上面tad和染色体热图填充一致
filenum<-332
filesize<-as.numeric(fs[filenum,1])
chrom<-'chr19'
data1<-0
data1<-read.table(filePath[filenum],header = F)
colnames(data1)<-c('chr1','pos1','chr2','pos2')
data1$chr1<-as.character(data1$chr1)
data1$chr2<-as.character(data1$chr2)
data1$interval1<-round(data1$pos1/resolution)
data1$interval2<-round(data1$pos2/resolution)
s1=0
s2=0
a1=0
a2=0
s1<-as.numeric(matsize_com[which(matsize_com[,1]==chrom),2])
s2<-as.numeric(matsize_com[which(matsize_com[,1]==chrom),3])
y<-matrix(0,nrow =s2,ncol =s2)
a1<-which(data1$chr1==chrom)
a2<-which(data1$chr2[a1]!=chrom)
if(length(a2)>0){
  data1<-data1[-a1[a2],]
}
if(length(which(data1$interval1[which(data1$chr1==chrom)]>s2))>0){
  data1<-data1[-which(data1$chr1==chrom)[which(data1$interval1[which(data1$chr1==chrom)]>s2)],]
}
if(length(which(data1$interval2[which(data1$chr1==chrom)]>s2))>0){
  data1<-data1[-which(data1$chr1==chrom)[which(data1$interval2[which(data1$chr1==chrom)]>s2)],]
}
for (i in which(data1$chr1==chrom)) 
{
  y[data1[i,5],data1[i,6]]<-1+y[data1[i,5],data1[i,6]];
}
for (i in which(data1$chr1==chrom)) 
{
  y[data1[i,6],data1[i,5]]<-1+y[data1[i,6],data1[i,5]];
}
y<-log((y/filesize)+1)
#y[which(y>0)[which.min(y[y>0])]]
#y[which.max(y)]

#下面这段代码用于计算：在基准文件中，该条染色体上的每个compartment有多少行，以用于处理可区分A/B compartment的热图矩阵
a1<-which(compartments$chrom==chrom)[1]#a1和a2记录了该条染色体在compartment基准文件里面的起始行和终止行
a2<-which(compartments$chrom==chrom)[length(which(compartments$chrom==chrom))]
i=0
j=1
n=0#该指针用于更新记录当前compartment的行数
num_of_com<-array()#该数组用于记录该条染色体上的每个compartment有多少行
for (i in a1:a2) {#该循环为计算compartment行数的方法
  if(i==a2){#若循环到最后一行，则表明无论如何都是一个compartment，将n填入对应的num_of_com并结束循环
    n=n+1
    num_of_com[j]<-n
    break()
  }
  if(compartments$Evalue[i]*compartments$Evalue[i+1]<0){#如果相邻两行的compartment值异号，则意味着compartment的切换，将n填入对应的num_of_com并清零
    n=n+1
    num_of_com[j]<-n
    j=j+1
    n=0
  }
  else{#如果相邻两行的compartment值同号，则继续累加
    n=n+1
  }
}

#我们在第一步就已经将compartment的热图矩阵绘制完成了，下面这段代码是利用上面得到的num_of_com，将热图矩阵当中的A/B compartment对应的数值用正负区分开来，这样在绘制热图的时候就能形成鲜明对比
compar<-compartments[which(compartments$chrom==chrom),]#将基准文件中，该条染色体对应的行取出来存到compar中
ncom<-0#这个指针指向各个compartment的行数：num_of_com[ncom]
ncp<-0#这个指针指向上一个compartment的最后一行：ncp=ncp+num_of_com[ncom]
for(ncom in 1:length(num_of_com)){#下面为将热图矩阵中的A/B compartment区分为正负数值的方法
  if(ncom==length(num_of_com)){#如果ncom指向最后一个compartment
    for (i in compar$start[ncp+1]:compar$end[ncp+num_of_com[ncom]]) {#将最后一个compartment对应的热图区间乘以对应的compartment值（+1/-1）
      y[i,compar$start[ncp+1]:ncol(y)]<-y[i,compar$start[ncp+1]:ncol(y)]*compar$Evalue[ncp+1]
    }
  }
  else{#如果ncom不指向最后一个compartment，则下面有三步操作
    for (i in compar$start[ncp+1]:(compar$end[ncp+num_of_com[ncom]]-1)) {
      y[i,(compar$end[ncp+num_of_com[ncom]]+1):ncol(y)]<-0#首先第一步操作是将不属于该compartment区间的行与列之外的数据删除
      y[(compar$end[ncp+num_of_com[ncom]]+1):ncol(y),i]<-0
    }
    for (i in compar$start[ncp+1]:compar$end[ncp+num_of_com[ncom]]) {#第二步是将该compartment对应的热图区间乘以对应的compartment值（+1/-1）
      y[i,compar$start[ncp+1]:compar$end[ncp+num_of_com[ncom]]]<-y[i,compar$start[ncp+1]:compar$end[ncp+num_of_com[ncom]]]*compar$Evalue[ncp+1]
    }
    #第三步是将该compartment与下一个compartment共交的那个矩阵位点归于下一个compartment，所以这一步需要将这一矩阵位点的数值再次乘以当前compartment对应的compartment值（+1/-1）以达到还原效果
    y[compar$end[ncp+num_of_com[ncom]],compar$end[ncp+num_of_com[ncom]]]<-y[compar$end[ncp+num_of_com[ncom]],compar$end[ncp+num_of_com[ncom]]]*compar$Evalue[ncp+1]
    ncp=ncp+num_of_com[ncom]#该层循环结束，指针更新
  }
}

#绘制可区分A/B compartment的热图
ramp <- colorRamp(c("firebrick3", "white","navy"));#这里设置颜色梯度为暖色-白色-冷色
cols<-rgb( ramp(seq(0, 1, length = 128)), max = 255);
tit<-paste('sample',filenum,'_',chrom,'(',startsite,'-',endsite,'/',resolution/1000,'k)',sep = '')
heatmap(y[compar$start[1]:ncol(y),compar$start[1]:ncol(y)],col=rev(cols), Rowv=NA,Colv=NA, labRow=NA,labCol=NA,scale="none",main=tit);#绘制热图的区间需要参照基准文件的尺度标准（matsize_com）
#下面为选择一小段区间的热图，注：d3heatmap为交互式热图，需要耗费大量资源，只能在小范围使用
ramp <- colorRamp(c("firebrick3", "white","navy"));
cols<-rgb( ramp(seq(0, 1, length = 128)), max = 255);
heatmap(y[87:185,87:185],col=rev(cols), Rowv=NA,Colv=NA, labRow=NA,labCol=NA,scale="none",main=tit);
d3heatmap(y[87:185,87:185], colors = "RdBu",Rowv = F,Colv = F,show_grid = T)












#chromosome counts of mc (之前尝试的处理甲基化的数据代码，效果较差，不用管它=_=#)
path<-"/store/yhshao/GSE124391_RAW/all_txt";
fileNames <- dir(path) ;
sana_mESC<-read.csv(file = '/store/yhshao/mESC.txt',header = T)
sana_mESC<-sana_mESC[order(sana_mESC$GEO_Accession..exp.),]
sana_mESC<-sana_mESC[-c(1:3),]
fileNames<-fileNames[-which(sana_mESC$Organism=='Homo sapiens')]
sana_mESC<-sana_mESC[-which(sana_mESC$Organism=='Homo sapiens'),]
qc<-array()
# j=1
# for (filenum in 1:length(filePath)) {
#   data1<-0
#   data1<-read.table(filePath[filenum],header = F)
#   colnames(data1)<-c('chr1','pos1','chr2','pos2')
#   data1$chr1<-as.character(data1$chr1)
#   data1$chr2<-as.character(data1$chr2)
#   qc[j]<-length(which(data1$chr1==data1$chr2))/nrow(data1)
#   j=j+1
#   print(paste0('第',filenum,'条样本统计完毕'))
# }
qc<-read.csv(file='/store/yhshao/Rwork/qc.csv')
qc<-qc[,-1]
fileNames<-fileNames[-which(qc<0.5)]
sana_mESC<-sana_mESC[-which(qc<0.5),]
filePath <- sapply(fileNames, function(x){ paste(path,x,sep='/')})
print(paste("共",length(filePath),"条待处理样本"))

t1=proc.time()
chrs<-seq(1,19)
chrs[20:21]<-c('X','Y')
samplematrix_mc_chr<-matrix(0,nrow = length(filePath),ncol = length(chrs))
for (filenum in 1:length(filePath)) {
  data1<-0
  data1<-read.table(filePath[filenum],header = F)
  chrsite<-1
  for (chrom in chrs) {
    samplematrix_mc_chr[filenum,chrsite]<-sum(data1$V5[which(data1$V5[which(data1$V1==chrom)]!=0)])
    chrsite=chrsite+1
  }
  print(paste('第',filenum,'/',length(filePath),'条样本计数完毕(mc_chr)'))
  
}
write.csv(samplematrix_mc_chr,file = '/store/yhshao/Rwork/samplematrix_mc_chr/samplematrix_mc_chr')
t2=proc.time()
t=t2-t1
print(paste0('分染色体甲基化计数执行时间：',t[3][[1]]/60,'分'))



