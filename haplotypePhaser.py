import numpy as np
import collections


#global constant
HAPLO_SINGLE_LIM = 512 # optimization
HAPLO_NUM_LIM = 2000 # optimization
WINDOW_SIZE_LIM = 18  # optimization, limit for the dynamic window size

OVERLAP_SNP_NUM = 6 # optimization, number of SNP overlapped between two windows
OVERLAP_SNP_DEC = 64 # 2^OVERLAP_SNP_NUM, e.g. 2^3 = 8

SEGMENT_START = 0 #从39074开始，节省时间，正常情况下调为0

#global variable
genoList=[] # (50*num_of_SNP, e.g. 50*39496)
indivQueue = collections.deque() # deque is "double-ended queue"
maxList=[] # list of the best haplotype pairs and corresponding p-value for each individual (50*3)
oldMaxList=[] # stores the last OVERLAP_SNP_NUM of last maxList (50*3)
newList=[] # oriented maxList after orient() (50*3)

#haplotype from segStart,itemNum
segStart=0;     # handle with a segment (window) starting with segStart,itemNum from genoList     
itemNum=WINDOW_SIZE_LIM;  

haplotypList = []  

idx_hapQueue = collections.deque()
idx_indivQueue = collections.deque()

pvalueInitList = []
pvalueNextList = []


#1: 读取原始数据，存储为genoList[]
def readGenoFile():
    global file_output

    print("Please type your file paths to start phasing.")
    file_input = input("Input file path (press Enter when finish): ")
    file_output = input("Output file name (press Enter when finish): ")
    # file_path=r'/Users/xinhuang/Documents/Junior/2019 Winter/CS CM 124/Programming Project/Assignment Details/test_data_1.txt'
    with open(file_input, "r") as file_object:
        lines=file_object.readlines() # read each line
    
    # print(lines[0]) # lines[0] is the first SNP for all 50 individuals

    nums=lines[0].split() # split to each individual

    for i in range(0,len(nums)): # create a matrix genoList[] and fill with '0'
        tempList=[]
        for g in range(0,len(lines)):
            tempList.append(0)
        genoList.append(tempList)

    #print(genoList[1][0])
    #print(genoList[0][1])
    
    row=0;
    col=0;
    for line in lines:
        nums=line.split()
        col=0;
        for num in nums:
            genoList[col][row]=int(num)
            col=col+1
        row=row+1
    #print('row='+str(row))
    #print('col='+str(col))
    nRows=row
    nCols=col

    num0=0;
    num1=0;
    num2=0;
    for i in range(0,row): # count the number of 0, 1, 2 in the first row ???
        if genoList[0][i]==1:
            num1=num1+1
        if genoList[0][i]==0:
            num0=num0+1
        if genoList[0][i]==2:
            num2=num2+1
    
    #print('num0='+str(num0))
    #print('num1='+str(num1))
    #print('num2='+str(num2))
    return
#end readGenoFile()

readGenoFile()

#2: 对于个体i，拆解一段genoList[]为所有可能的haplotypes，存储为haplotypList[]
def addlist(hapList,individual_i,segStart,itemNum):
#    individual_i=1;  #individual
# segStart is the start of the window (e.g. 1-32, 4-36)
# itemNum is the window size

    global genoList
    global indivQueue

    tempList = [] 
    tempList.append(0)
    #for loop to establish all the haplotypes
    for j in range(0,itemNum):
        if genoList[individual_i][j+segStart]==2:
            for k in range(0,len(tempList)):
                tempList[k]=tempList[k]+(1<<j)  
                # little-ending(小端模式) 处理number时，从最低位开始(右)，比如‘321’首先处理‘1’，而不是‘3’；和字符串string是完全反过来的
                # tempList[k] is a number, instead of a array, 不是一般的位操作
                # i<<j means 左移j位, 低位在右，高位在左
        elif genoList[individual_i][j+segStart]==1:
            tempList.extend(tempList[:]) 
            # tempList[:]表示全集，extend(tempList[:])就是扩充一倍
            # array.extend(array[n:m])扩充array从n到m的所有entry附加在array之后
#            print(len(tempList))
            if len(tempList)>HAPLO_SINGLE_LIM:
                return -1
            for k in range(len(tempList)//2,len(tempList)): #Python里‘//’是整数除，‘/’是浮点数除
                tempList[k]=tempList[k]+(1<<j)
#    print(len(tempList))
    indivQueue.append(tempList[:])

    #合并“分支hap”列表到“主hap”列表中
    hapList.extend(tempList) # append是加一个元素，extend是加一个集合
    
    #排序，以利查重
    haplotypList.sort()
    
    #查找到重复串，删除
    for j in range(0,len(hapList)-1): # In Python, 自动减一了，这里要额外减一个1
        if hapList[j]==hapList[j+1]:        
            del hapList[j+1]
        if j>=len(hapList)-2:
            break

    return len(haplotypList)
#end addlist()


#3: 根据长短决定window size，并挨个individual执行addlist()
def findWinSize():
    global segStart,itemNum
    global indivQueue
    global haplotypList

    low=0
    high=WINDOW_SIZE_LIM
    while True:
        indivQueue.clear()
        haplotypList.clear()
        
        itemNum=min((low+high)//2,len(genoList[0])-segStart)
        for i in range(0,len(genoList)):    #50 for example1   
            hapLen=addlist(haplotypList,i,segStart,itemNum)
            if hapLen==-1:
                break
        if (hapLen==-1) or (hapLen>HAPLO_NUM_LIM) :
            if high==itemNum:
                break
            high=itemNum
        else:
            if low==itemNum:
                break
            low=itemNum      
    
    # print("segStart:"+str(segStart))
    # print("hapLen:"+str(hapLen)+";last itemNum:"+str(itemNum)+";low:"+str(low)+";high:"+str(high))

    return 
#end findWinSize()


#4: 计算对偶
def mate(indiv,mateone):
    global genoList
    global segStart,itemNum

    mateother=0
    for j in range(0,itemNum):
        if genoList[indiv][j+segStart]==2:
            mateother=mateother+(1<<j)
        elif genoList[indiv][j+segStart]==1:
            if (mateone &(1<<j))==0: 
                mateother=mateother+(1<<j)
    return mateother
#end mate()


#5: 建立hap索引
def buid_idx_hapQueue():
    global haplotypList
    global indivQueue
    global idx_indivQueue
    global idx_hapQueue

    idx_hapQueue.clear()

    for k in range(0,len(haplotypList)):
        tempList = [] 
        for i in range(0,len(indivQueue)):
            try: # try lets you test a block of code for errors，因为indivQueue[i]里面不一定有haplotypList[k]
                idx=indivQueue[i].index(haplotypList[k])
                # 交叉索引，indivQueue[i]直接就是类对象，在indivQueue[i]里面找到haplotypList[k]的索引
            except: # when the try error occurs, run the except，意外应该经常发生
                continue 
            tempList.append(i) # 他自己就是自己的索引
            idxm=haplotypList.index(mate(i,haplotypList[k]))   #must exist!
            tempList.append(idxm) # 存储对偶的那个索引
        idx_hapQueue.append(tempList[:])
#end buid_idx_hapQueue()


#6: 建立indivQueue索引
def build_idx_indivQueue():
    global haplotypList
    global indivQueue
    global idx_indivQueue

    idx_indivQueue.clear()

    for i in range(0,len(indivQueue)):
        tempList = [] 
        #print('ia'+str(i)+'='+str(len(indivQueue[i])))
        for j in range(0,len(indivQueue[i])):
            if j>len(indivQueue[i])-1:
                break
            kkk=haplotypList.index(indivQueue[i][j])
            xxx=mate(i,indivQueue[i][j])
            lll=haplotypList.index(xxx)
            tempList.append(kkk)
            tempList.append(lll)
            if kkk!=lll:
                indivQueue[i].remove(xxx)
                #indivQueue中只保留对偶（k,l)中，数较小的，也是list里靠前的，删掉后者，但要把二者索引都找到，并记下；
        
        idx_indivQueue.append(tempList)
    
#end build_idx_indivQueue()


#7: 建立？？？
def build_pvalue_indivQueue():
    global pvalue_indivQueue
    global indivQueue
    global maxList

    maxList.clear()
    pvalue_indivQueue.clear()

    for i in range(0,len(indivQueue)):
        tempList = [] 
        maxItem=[0,1,0.9]   #left haplo(0),right hiplo(1) of maxItem,3 maxPvalue(0.9),init any value,OK!
        for j in range(0,len(indivQueue[i])):
            tempList.append(0)
        pvalue_indivQueue.append(tempList)
        maxList.append(maxItem)
#end build_pvalue_indivQueue()


#8: 初始化pvalue为平均分配：1/(len(haplotypList))，顺序和haplotypList一致
def build_pvalueInitList():
    for h in range(0,len(haplotypList)):
        pvalueInitList.append(1/(len(haplotypList)))
#end build_pvalueInitList()


#9: 计算更新pvalue，算法的精华
def updateP(pvalueList):
#calculate denominator,sum(pk*pl) for indiv i
#    print('denomList')
    global haplotypList
    global idx_indivQueue
    global idx_hapQueue
    global pvalue_indivQueue

    denomList = []
    for i in range(0,len(idx_indivQueue)):
        denom=0
        for kl in range(0,len(idx_indivQueue[i])//2):
            kkk=idx_indivQueue[i][kl*2]     #kkk,index to left one of mate(pair)
            lll=idx_indivQueue[i][kl*2+1]   #kkk,index to right one of mate(pair)
            denom=denom+pvalueList[kkk]*pvalueList[lll]                
        denomList.append(denom)
#        print(str(denom))
    
    #new value calc
    pvalueNewList=[]
    nnn=len(indivQueue)
    for h in range(0,len(haplotypList)):
        newValue=0
        for ik in range(0,len(idx_hapQueue[h])//2):
            iii=idx_hapQueue[h][ik*2]       #iii,index to individual,(from 0 to 49,in example_data_1)
            kkk=idx_hapQueue[h][ik*2+1]     #kkk,index to mate of h
            if h==kkk:
                newValue=newValue+2*pvalueList[h]*pvalueList[kkk]/denomList[iii]
            else:
                newValue=newValue+pvalueList[h]*pvalueList[kkk]/denomList[iii]
                
        pvalueNewList.append(newValue/(2*nnn))
    
#        print(str(pvalueNewList[h]))

    #using pvalueNewList,to compute probability of mate
    for i in range(0,len(idx_indivQueue)):
        sum=0;
        for kl in range(0,len(idx_indivQueue[i])//2):
            kkk=idx_indivQueue[i][kl*2]     #kkk,index to left one of mate(pair)
            lll=idx_indivQueue[i][kl*2+1]   #kkk,index to right one of mate(pair)
            sum=sum+pvalueNewList[kkk]*pvalueNewList[lll]
            pvalue_indivQueue[i][kl]=pvalueNewList[kkk]*pvalueNewList[lll]
            
        for kl in range(0,len(idx_indivQueue[i])//2):
            pvalue_indivQueue[i][kl]=pvalue_indivQueue[i][kl]/sum
        
        kl=np.argmax(pvalue_indivQueue[i])
        kkk=idx_indivQueue[i][kl*2]
        lll=idx_indivQueue[i][kl*2+1]
        maxList[i][0]=haplotypList[kkk]
        maxList[i][1]=haplotypList[lll]
        maxList[i][2]=pvalue_indivQueue[i][kl]
        
        pvalue_indivQueue[i].sort(reverse=True)

    #test...    
    # sumValue=0
    # for h in range(0,len(haplotypList)):
    #     sumValue=sumValue+pvalueNewList[h]
    # print('sum of pvalueNew (should be 1.0)='+str(sumValue))
    return pvalueNewList
#end updateP()


#10: 最后关闭文件（必须关闭，关闭才能写到硬盘上）
def fileEnd():
    global file_object
    file_object.close()
#end fileEnd()


#11: 决定正反顺序
def orient():
    global maxList
    global oldMaxList
    global newList

    # newList.clear()
    # print(newList[5][1])

    for i in range(0,len(idx_indivQueue)):
        newList[i][0]=maxList[i][0]
        newList[i][1]=maxList[i][1] # 初始化

    if (segStart-SEGMENT_START > OVERLAP_SNP_NUM): # 如果不是第一个segment
        for i in range(0,len(idx_indivQueue)):
            if (maxList[i][0] % OVERLAP_SNP_DEC)==oldMaxList[i][1]:
                newList[i][0]=maxList[i][1]
                newList[i][1]=maxList[i][0] # 反着对得上，前后调换
            #     print("reverse")
            # elif (maxList[i][0] % OVERLAP_SNP_DEC)!=oldMaxList[i][0]:
            #     if pvalue_indivQueue[i][0]-pvalue_indivQueue[i][1] < 0.15:
            #         # 取前一序列的后三个
            #     else:
            #         # 取当前序列的前三个
            #     print(segStart, "BADindiv", i, pvalue_indivQueue[i][0:2])
                # print("first 3 maxList", maxList[i][0] % OVERLAP_SNP_DEC)
                # print("last 3 oldMaxList", oldMaxList[i][0])
            # 本来就对得上 或者 完全对不上，都原地不动
            newList[i][2]=maxList[i][2]
#end orient()


#12: 用newList内容输出文件
def fileOut():
    global file_object

    # for i in range(0,len(idx_indivQueue)):
    #     print(pvalue_indivQueue[i][0:3]) #打印最大的三个pvalue
    
    rowStr=[]

    if (segStart-SEGMENT_START==0 or itemNum<=OVERLAP_SNP_NUM): # 如果是第一个或者最后一个segment
        for h in range(0,itemNum):  #itemNum is haplo length in this case
            s0=''
            for i in range(0,len(idx_indivQueue)):
                if (newList[i][0]&(1<<h))==0:
                    s0=s0+"0 "  #add space
                else: 
                    s0=s0+"1 "
                if (newList[i][1]&(1<<h))==0:
                    s0=s0+"0"   #NO space
                else: 
                    s0=s0+"1"
                if i==(len(idx_indivQueue)-1):
                    s0=s0+"\n"
                else:
                    s0=s0+" "
        
            rowStr.append(s0)
            # print(rowStr[h],end='')

        for h in range(0,itemNum):  #itemNum is haplo length for this case
            file_object.write(rowStr[h]) # 写入

    else: # 如果不是第一个segment
        for h in range(OVERLAP_SNP_NUM,itemNum):  #itemNum is haplo length in this case
            s0=''
            for i in range(0,len(idx_indivQueue)):
                if (newList[i][0]&(1<<h))==0:
                    s0=s0+"0 "  #add space
                else: 
                    s0=s0+"1 "
                if (newList[i][1]&(1<<h))==0:
                    s0=s0+"0"   #NO space
                else: 
                    s0=s0+"1"
                if i==(len(idx_indivQueue)-1):
                    s0=s0+"\n"
                else:
                    s0=s0+" "
        
            rowStr.append(s0)
            # print(rowStr[h-OVERLAP_SNP_NUM],end='')
    
        for h in range(OVERLAP_SNP_NUM,itemNum):  #itemNum is haplo length for this case
            file_object.write(rowStr[h-OVERLAP_SNP_NUM]) # 写入
#end fileOut()


#13: 储存当次maxList的最后OVERLAP_SNP_NUM位数字，以备下一循环调用
def storeList():
    global newList
    global oldMaxList

    # oldMaxList.clear()
    for i in range(0,len(idx_indivQueue)):
        oldMaxList[i][0]=newList[i][0]>>(itemNum-OVERLAP_SNP_NUM) # 取二进制数最高OVERLAP_SNP_NUM位
        # print("MaxList_num: ", oldMaxList[i][0])
        oldMaxList[i][1]=newList[i][1]>>(itemNum-OVERLAP_SNP_NUM)
        oldMaxList[i][2]=newList[i][2]
#end storeList()

########################################################################################
#Main Program!
#global variable clr

#haplotype from segStart,itemNum   
segStart=SEGMENT_START;     #handle with a segment starting with segStart,itemNum from genoList     
itemNum=WINDOW_SIZE_LIM;  

pvalueInitList = []
pvalueNextList = []

pvalue_indivQueue = collections.deque()
maxList=[]
oldMaxList = [[0 for col in range(3)] for row in range(len(genoList))] # stores the last OVERLAP_SNP_NUM of last maxList
newList = [[0 for col in range(3)] for row in range(len(genoList))]

import time;  # 引入time模块
oldticks = time.time()
print("Current time is:", oldticks)

# fileout_path=r'test_data_1_sol.txt'
fileout_path=file_output
file_object=open(fileout_path,'w')
file_object.close()
file_object=open(fileout_path,'a+') 

while True: #循环，直到文件末尾
    indivQueue.clear
    idx_hapQueue.clear
    idx_indivQueue.clear
    pvalue_indivQueue.clear()
    maxList.clear()
    pvalueInitList.clear()

    findWinSize()
    ticks = time.time()

    # idx_hapQueue building
    # print('hap='+str(len(haplotypList)))
    # print('indiv='+str(len(indivQueue)))

    buid_idx_hapQueue()
    build_idx_indivQueue()
    build_pvalue_indivQueue()
    build_pvalueInitList()

    # EM iteration, 迭代三次就能收敛
    pvalueNextList=updateP(pvalueInitList)
    pvalueNextList=updateP(pvalueNextList)
    pvalueNextList=updateP(pvalueNextList)

    orient() # switch haplo according to overlap
    fileOut()
    storeList() # store the old haplo list

    if itemNum>OVERLAP_SNP_NUM:
        segStart=segStart+itemNum-OVERLAP_SNP_NUM # overlap!!!
    else: #itemNum<=OVERLAP_SNP_NUM, near end of file
        segStart=segStart+itemNum

    print("Running Time:", round(ticks-oldticks, 1), "s")
    print("Speed:       ", round((segStart-SEGMENT_START)/(ticks-oldticks), 1), "SNP/s")
    print("Progress:    ", 100*(segStart-SEGMENT_START)//len(genoList[0]), "%")

    if segStart>=(len(genoList[0])-OVERLAP_SNP_NUM): # end of file
        break

fileEnd()

print("OK end!!!")
input()



