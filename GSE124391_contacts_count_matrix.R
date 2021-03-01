#初始化文件信息
path<-"/store/yhshao/GSE124391_RAW";#交互信息文件的路径
resolution<-100000#全局分辨率
fileNames<-dir(path)#按照路径获取文件名
fileNames<-fileNames[-c(1:5)]#去除前五个（多细胞Hi-C数据以及其他不相关文件）
sana_mESC<-read.csv(file = '/store/yhshao/mESC.txt',header = T)#读入样本详细信息文件
sana_mESC<-sana_mESC[order(sana_mESC$GEO_Accession..exp.),]#对样本信息文件按照GEO号进行排序，因为路径里面的样本是按照GEO号排序的
sana_mESC<-sana_mESC[-c(1:3),]#去除样本信息文件中多细胞Hi-C样本的信息
fileNames<-fileNames[-which(sana_mESC$Organism=='Homo sapiens')]#将样本列表和样本信息文件中的人类样本去除
sana_mESC<-sana_mESC[-which(sana_mESC$Organism=='Homo sapiens'),]


#提取细胞类型名称
qc<-array()#qc是一个数组，用于存储每个样本的染色体内交互占比，正常情况下一个单细胞数据的染色体内交互占比应该为50～60%
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

#下面为染色体、合并与非合并A/B compartments三个层面的计数矩阵算法
t1=proc.time()#记录程序运行开始时间
samplematrix_chr<-matrix(0,nrow = length(filePath),ncol=length(chrs))#染色体层面的计数矩阵，行代表样本，列代表染色体：358*21
samplematrix_com<-matrix(0,nrow = length(filePath),ncol=length(chrs)*2)#合并compartment层面的计数矩阵，行代表样本，列代表染色体的A/B两个compartment：358*42
samplematrix_com2<-matrix(0,nrow = length(filePath),ncol=sum(as.numeric(matsize_com[,4])))#非合并compartment层面的计数矩阵，行代表样本，列代表单独的一个compartment：358*n(compartment)
for (filenum in 1:length(filePath)) {#最外层大循环，每次处理一个样本文件
  data1<-0#初始化文件存储变量
  data1<-read.table(filePath[filenum],header = F)#读入文件为表格
  colnames(data1)<-c('chr1','pos1','chr2','pos2')#将列名称填入
  data1$chr1<-as.character(data1$chr1)#将1、2号染色体名称的数据型改为字符串
  data1$chr2<-as.character(data1$chr2)
  data1$interval1<-round(data1$pos1/resolution)#将交互位点的数值除以分辨率并四舍五入，存入interval列
  data1$interval2<-round(data1$pos2/resolution)
  a1<-0
  a2<-0
  chrsite<-1#在染色体层面算法中指向计数矩阵里面每一列（染色体）的指针
  intra_contacts<-0#在染色体层面算法中记录染色体内交互的行信息
  numchrom<-1#在合并compartment算法中指向计数矩阵里面每一列（每条染色体的A/B compartment）的指针
  comparsite<-1#在非合并compartment算法中指向计数矩阵里面每一列（每个单独compartment）的指针
  
  for (j in chrs) {#这一个for循环用于记录染色体层面的交互计数
    a1<-which(data1$chr1==j)#a1和a2用于记录1、2号染色体名称一致的行
    a2<-which(data1$chr2[a1]==j)
    intra_contacts<-a1[a2]#将该染色体内交互的行记录下来
    if(length(intra_contacts)!=0){#如果该染色体的染色体内交互行数不为零，则将行数记为该条染色体在染色体层面的交互计数
      samplematrix_chr[filenum,chrsite]=samplematrix_chr[filenum,chrsite]+length(intra_contacts)
      chrsite=chrsite+1
    }
    else{#如果该染色体的染色体内交互行数为零，则证明该文件没有该染色体的染色体内交互数据，将指向计数矩阵染色体的指针移向下一位
      chrsite=chrsite+1
    }
  }
  
  for (chrom in chrs) {#这一个for循环用于记录合并与非合并compartment层面的交互计数
    if(length(which(data1$chr1==chrom))!=0){#首先判断文件是否有该条染色体交互信息，如果有则继续操作
      s1=0
      s2=0
      a1=0
      a2=0
      compar<-0
      compar<-compartments[which(compartments$chrom==chrom),]#将基准compartment文件中属于该染色体的行取出来存入compar
      compar$count<-compar$Evalue*0#在compar中生成一个新列count，用于存储每个compartment的交互计数
      s1<-as.numeric(matsize_com[which(matsize_com[,1]==chrom),2])#s1和s2用于存储绘制该染色体热图的矩阵大小（注：热图的绘制尺度是基于基准文件的（matsize_com））
      s2<-as.numeric(matsize_com[which(matsize_com[,1]==chrom),3])
      y<-matrix(0,nrow =s2,ncol =s2)#初始化热图矩阵
      a1<-which(data1$chr1==chrom)#找出该样本文件中第一号染色体名称对应当前染色体名称的行
      a2<-which(data1$chr2[a1]!=chrom)#在上述行（a1）中找出第二号染色体名称不等于第一号染色体名称的行
      if(length(a2)>0){
        data1<-data1[-a1[a2],]#若a2不为零，那么代表存在染色体间交互数据，将其去除
      }
      if(length(which(data1$interval1[which(data1$chr1==chrom)]>s2))>0){#找到样本文件中第一个交互位点（interval1）的尺度超过基准文件尺度的行并去除
        data1<-data1[-which(data1$chr1==chrom)[which(data1$interval1[which(data1$chr1==chrom)]>s2)],]
      }
      if(length(which(data1$interval2[which(data1$chr1==chrom)]>s2))>0){#找到样本文件中第二个交互位点（interval2）的尺度超过基准文件尺度的行并去除
        data1<-data1[-which(data1$chr1==chrom)[which(data1$interval2[which(data1$chr1==chrom)]>s2)],]
      }
      for (i in which(data1$chr1==chrom)) #使用文件中剩余的染色体内交互填充热图矩阵
      {
        y[data1[i,5],data1[i,6]]<-1+y[data1[i,5],data1[i,6]];
      }
      for (i in which(data1$chr1==chrom)) 
      {
        y[data1[i,6],data1[i,5]]<-1+y[data1[i,6],data1[i,5]];
      }
      
      i=1#该指针指向compar的每一行
      j=0#该变量用于记录compartment变换时所经历的行数
      s=compar$end[i]#该变量记录第i行的终止位点
      counts=0#该变量用于记录每个单独compartment区域内的交互计数
      while (i!=nrow(compar)+1) {#对compar进行遍历
        if(i==nrow(compar)){#如果i指针指向最后一行
          s=compar$end[i]#更新s变量
          for (m in compar$end[i]:compar$start[(i-j)]) {#这里的双重for循环是为了对热图矩阵里面这个compartment区域的交互进行计数，
            for (n in compar$start[(i-j)]:s) {          #由于热图是对称的，我们需要计数的只有对角线和对角线一侧的数值，这里的m变量是
              counts=counts+y[m,n]                      #从第i行的终止位点到这个compartment的第一行的起始位点，n变量是从这个compartment的
            }                                           #第一行的起始位点到s变量并且内层循环结束之后s变量减1，所以，假设这个compartment的
            s=s-1                                       #区间为1～3，那么通过m与n就会生成（3，1）（3，2）（3，3）（2，1）（2，2）（1，1）的序列
          }
          compar$count[i-j]<-counts#将这个compartment的交互计数存入这个compartment的第一行（i-j）的count列
          samplematrix_com2[filenum,comparsite]<-counts#将这个compartment的交互计数存入非合并compartment计数矩阵的相应位置
          comparsite=comparsite+1#指向下一个compartment
          counts=0#计数归零
          i=i+1#指向compar的下一行
          j=0#当前compartment行数计数归零
          break()#对最后一个compartment计数完成，直接退出循环
        }
        if(compar$Evalue[i]*compar$Evalue[i+1]<0){#若i指针不指向最后一行，且当前i指针指向的行，其compartment值与相邻的下一行compartment值异号，则说明发生了A/B compartment的切换
          s=compar$end[i]
          for (m in compar$end[i]:compar$start[(i-j)]) {#这里的compartment区域计数与上述相同
            for (n in compar$start[(i-j)]:s) {
              counts=counts+y[m,n]
            }
            s=s-1
          }
          #由于这里不是最后一个compartment，因此会存在compartment交界处位点归属的问题，例如：一个A compartment的区间为1～3，与其相邻的下一个B compartment
          #的区间为3～4，那么3就为两者所共有的位点，在这里我将共有位点统一归为相邻两者的后者所拥有，也就是说在热图矩阵上（3，3）这个位点认为属于B compartment，
          #但是上面的双重循环还是按照1～3，3～4区间进行计数的，因此我们在对每一个非末尾compartment（A compartment的1～3）的区间进行双重循环计数之后，都要减去
          #交界处位点（3，3）的数值
          compar$count[i-j]<-counts-y[compar$end[i],compar$end[i]]#减去交界位点数值并存入这个compartment的第一行（i-j）的count列
          samplematrix_com2[filenum,comparsite]<-counts-y[compar$end[i],compar$end[i]]#将这个compartment的交互计数减去交界位点数值并存存入非合并compartment计数矩阵的相应位置
          comparsite=comparsite+1                                                     
          counts=0
          i=i+1
          j=0
        }
        else{#既不是compartment交换点也不是最后一个compartment的情况
          i=i+1
          j=j+1
        }
      }
      #while循环结束，现在只需要将compar中count列进行加和，就可以得到合并compartment的结果了
      samplematrix_com[filenum,numchrom]<-sum(compar$count[which(compar$Evalue>0)])#将compartment值大于零（A compartment）的count加和并存入合并compartment计数矩阵的对应位置
      samplematrix_com[filenum,numchrom+1]<-sum(compar$count[which(compar$Evalue<0)])#将compartment值小于零（B compartment）的count加和并存入合并compartment计数矩阵的对应位置
      numchrom=numchrom+2
    }
    else{#承接大循环里面第一个if语句，如果文件没有该染色体的交互信息
      numchrom=numchrom+2#将指向计数矩阵染色体的指针加2（一条染色体合并为A/B compartment）
      comparsite=comparsite+as.numeric(matsize_com[which(matsize_com[,1]==chrom),4])-1#将指向计数矩阵单独compartment的指针跳转到下一条染色体的第一个compartment
    }
  }
  print(paste('第',filenum,'/',length(filePath),'条样本计数完毕(chr mergecom notmergecom)'))#每层大循环结束之后的提示
}

#上面的算法得到了染色体、合并/非合并compartment三个层面的计数矩阵，接下来是对得到的三个计数矩阵进行去除零列（全部为零的列无法进行降维）和存储操作
#对染色体层面计数矩阵的操作
zerocol<-0#该数组用于存储计数矩阵里面全部为零的列
j=1
for (i in 1:ncol(samplematrix_chr)) {#对计数矩阵每一列进行遍历
  if(sum(samplematrix_chr[,i])==0){
    zerocol[j]<-i
    j=j+1
  }
}
if(sum(zerocol)>0){#如果存在零列，将其报告并去除
  print('染色体层面计数矩阵中的零列为：')
  print(zerocol)
  samplematrix_chr<-samplematrix_chr[,-zerocol]
}else{print('染色体层面计数矩阵中没有零列')}
write.csv(samplematrix_chr,file = '/store/yhshao/Rwork/samplematrix_100k_ES_mESC_chr')#存储计数矩阵

#合并compartment层面
zerocol2<-0
j=1
for (i in 1:ncol(samplematrix_com)) {
  if(sum(samplematrix_com[,i])==0){
    zerocol2[j]<-i
    j=j+1
  }
}
if(sum(zerocol2)>0){
  print('合并compartment层面计数矩阵中的零列为：')
  print(zerocol2)
  samplematrix_com<-samplematrix_com[,-zerocol2]
}else{print('合并compartment层面计数矩阵中没有零列')}
write.csv(samplematrix_com,file = '/store/yhshao/Rwork/samplematrix_merge_100k_ES_mESC_com')

#非合并compartment层面
zerocol3<-0
j=1
for (i in 1:ncol(samplematrix_com2)) {
  if(sum(samplematrix_com2[,i])==0){
    zerocol3[j]<-i
    j=j+1
  }
}
if(sum(zerocol3)>0){
  print('非合并compartment层面计数矩阵中的零列为：')
  print(zerocol3)
  samplematrix_com2<-samplematrix_com2[,-zerocol3]
}else{print('非合并compartment层面计数矩阵中没有零列')}
write.csv(samplematrix_com2,file = '/store/yhshao/Rwork/samplematrix_notmerge_100k_ES_mESC_com')
t2=proc.time()#记录程序运行结束时间并打印
t=t2-t1
print(paste0('分染色体、compartment(合并、非合并)执行时间：',t[3][[1]]/60,'分'))


#TAD层面的计数矩阵算法
t1=proc.time()#记录程序运行起始时间
samplematrix_tad<-matrix(0,nrow = length(filePath),ncol=nrow(tads)-nrow(matsize_tad))#初始化TAD计数矩阵，由于TAD基准文件为相邻两行记录一个TAD区间的格式，因此TAD总数为基准文件行数减去染色体数
#315～354与上面三个层次的处理方式一致，都是读入文件并进行热图矩阵填充操作，之所以分开运行是因为matsize_com与matsize_tad所规定的热图尺度有所不同，并且当分辨率很大的时候会出现没有TAD的情况
for (filenum in 1:length(filePath)) {
  data1<-0
  data1<-read.table(filePath[filenum],header = F)
  colnames(data1)<-c('chr1','pos1','chr2','pos2')
  data1$chr1<-as.character(data1$chr1)
  data1$chr2<-as.character(data1$chr2)
  data1$interval1<-round(data1$pos1/resolution)
  data1$interval2<-round(data1$pos2/resolution)
  tadsite=1#该指针用于指向计数矩阵里面每一列（每一个TAD）
  for (chrom in chrs) {
    if(as.numeric(matsize_tad[which(matsize_tad[,1]==chrom),2])>0){#判断基准文件在这条染色体上是否有TAD，如果有，则继续进行操作
      if(length(which(data1$chr1==chrom))!=0){#判断当前文件是否有该条染色体信息，如果有，则继续进行操作
        s1=0
        s2=0
        a1=0
        a2=0
        compar<-0
        compar<-tads[which(tads$chrom==chrom),]#将TAD基准文件中该染色体的行取出存入compar
        s1<-as.numeric(matsize_tad[which(matsize_tad[,1]==chrom),2])#s1和s2用于存储绘制该染色体热图的矩阵大小（注：热图的绘制尺度是基于基准文件(matsize_tad)的）
        s2<-as.numeric(matsize_tad[which(matsize_tad[,1]==chrom),3])
        y<-matrix(0,nrow =s2,ncol =s2)
        a1<-which(data1$chr1==chrom)#找出该样本文件中第一号染色体名称对应当前染色体名称的行
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
        
        counts=0#该变量用于记录每个单独TAD区域内的交互计数
        for (i in 1:(nrow(compar)-1)) {#从第一行至倒数第二行遍历compar，i指针指向每一行
          if(i==(nrow(compar)-1)){#TAD基准文件是相邻两行就代表一个TAD区间，如果指针指向倒数第二行，则证明只剩下最后一个TAD需要进行计数了
            if(i==1){#如果i=1，证明基准文件在该条染色体上只有一个TAD且其只有一行数据
              #对该TAD区间进行双重循环加和的方法与之前所述类似，唯一不同点是，TAD基准文件是相邻两行就代表一个TAD区间，因此第一个TAD的区间取值
              #应该为compar$start～compar$end，而其余的TAD的区间取值为compar$end～compar$end，具体直接查看compar较为直观
              s=compar$end[i+1]
              for (m in compar$end[i+1]:compar$start[i]) {
                for (n in compar$start[i]:s) {
                  counts=counts+y[m,n]
                }
                s=s-1
              }
              samplematrix_tad[filenum,tadsite]<-counts#将该TAD区域内的计数存入计数矩阵相应位置，因为只有一个TAD，因此没有去除交界处位点的操作
              tadsite=tadsite+1#指针移向下一个TAD
            }
            else{#如果指针指向倒数第二行（最后一个TAD）且其不等于1，那将按照正常计算即可
              s=compar$end[i+1]
              for (m in compar$end[i+1]:compar$end[i]) {
                for (n in compar$end[i]:s) {
                  counts=counts+y[m,n]
                }
                s=s-1
              }
              samplematrix_tad[filenum,tadsite]<-counts#同理，最后一个TAD没有去除交界处位点的操作
              tadsite=tadsite+1
            }
            break()#基准文件只有一个TAD，计数完之后就可以结束循环了
          }
          if(i==1){#能够来到这一步，那就证明目前基准文件里至少有两个TAD，并且当前指针i指向第一行（第一个TAD）
            s=compar$end[i+1]
            for (m in compar$end[i+1]:compar$start[i]) {
              for (n in compar$start[i]:s) {
                counts=counts+y[m,n]
              }
              s=s-1
            }
            samplematrix_tad[filenum,tadsite]<-counts-y[compar$end[i+1],compar$end[i+1]]#因为目前基准文件里至少有两个TAD，需要进行交界位点去除的操作
            tadsite=tadsite+1
            counts=0
          }
          else{#指针不指向倒数第二行（最后一个TAD），也不指向第一行，那就证明目前基准文件里至少有三个TAD，并且目前指针i指向非两端的TAD中
            s=compar$end[i+1]
            for (m in compar$end[i+1]:compar$end[i]) {
              for (n in compar$end[i]:s) {
                counts=counts+y[m,n]
              }
              s=s-1
            }
            samplematrix_tad[filenum,tadsite]<-counts-y[compar$end[i+1],compar$end[i+1]]#因为目前基准文件里至少有两个TAD，也需要进行交界位点去除的操作
            tadsite=tadsite+1
            counts=0
          }
        }
      }
    }
    else{#承接大循环里面第一个if语句，如果基准文件在该染色体上没有TAD信息
      tadsite=tadsite+as.numeric(matsize_tad[which(matsize_tad[,1]==chrom),4])-1#将指向计数矩阵每个TAD的指针跳转到下一条染色体的第一个TAD位置
    }
  }
  print(paste('第',filenum,'/',length(filePath),'条样本计数完毕(TADs)'))#每层循环结束后的提示
}

#同上，下面是对零列的去除和存储操作
zerocol<-0
j=1
for (i in 1:ncol(samplematrix_tad)) {
  if(sum(samplematrix_tad[,i])==0){
    zerocol[j]<-i
    j=j+1
  }
}
if(sum(zerocol)>0){
  print('TAD层面计数矩阵中的零列为：')
  print(zerocol)
  samplematrix_tad<-samplematrix_tad[,-zerocol]
}else{print('TAD层面计数矩阵中没有零列')}
write.csv(samplematrix_tad,file = '/store/yhshao/Rwork/samplematrix_100k_ES_mESC_tad')
t2=proc.time()
t=t2-t1
print(paste0('分TAD执行时间：',t[3][[1]]/60,'分'))


