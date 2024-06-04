#比较同一张衍射图在三个volume中得到的cc max,然后进行分类
import numpy as np
import h5py
import sys, getopt

infile_cc_1=""
infile_cc_2=""
infile_cc_3=""
cc_dif=0 #统一张衍射图在两个volume下ccmax差异达到cc_dif以上,则认为可以被划分到cc更高的volume中
outfile_cc_1=""
outfile_cc_2=""
outfile_cc_3=""

#命令行输入
opts, args = getopt.getopt(sys.argv[1:],"h",["infile_cc_1=", "infile_cc_2=","infile_cc_3=", "cc_dif=", "outfile_cc_1=", "outfile_cc_2=", "outfile_cc_3="])
for opt, arg in opts:
    if opt in ("--infile_cc_1"):
        infile_cc_1 = arg
    elif opt in ("--infile_cc_2"):
        infile_cc_2 = arg
    elif opt in ("--infile_cc_3"):
        infile_cc_3 = arg
    elif opt in ("--cc_dif"):
        cc_dif = float(arg)
    elif opt in ("--outfile_cc_1"):
        outfile_cc_1 = arg
    elif opt in ("--outfile_cc_2"):
        outfile_cc_2 = arg
    elif opt in ("--outfile_cc_3"):
        outfile_cc_3 = arg
print("Compare CC begins...\n")

#对比cc系数
h5_incc_1=h5py.File(infile_cc_1,"r")
h5_incc_2=h5py.File(infile_cc_2,"r")
h5_incc_3=h5py.File(infile_cc_3,"r")

input_cc_1=h5_incc_1["angle"][:,:]  #上一轮所有衍射图与volume_1切片计算得到的cc结果
input_cc_2=h5_incc_2["angle"][:,:]
input_cc_3=h5_incc_3["angle"][:,:]

n_pattern=int(max(input_cc_1[:,0]))+1

output_cc_1=[] #所有被分类到第一个volume中的衍射图,下一轮用于merge volume_1
output_cc_2=[] 
output_cc_3=[] 

for i in range(n_pattern):
    cc_i_list=[input_cc_1[i,4],input_cc_2[i,4],input_cc_3[i,4]]
    cc_i_sorted = sorted([input_cc_1[i,4],input_cc_2[i,4],input_cc_3[i,4]])
    if cc_i_sorted[2]-cc_i_sorted[1]>cc_dif:
        max_cc=cc_i_sorted[2]
        index_max_cc=cc_i_list.index(cc_i_sorted[2])+1
        if index_max_cc==1:
            output_cc_1.append(input_cc_1[i,:])
        elif index_max_cc==2:
            output_cc_2.append(input_cc_2[i,:])
        elif index_max_cc==3:
            output_cc_3.append(input_cc_3[i,:])

#输出结果
h5_outcc_1=h5py.File(outfile_cc_1,"w")
h5_outcc_2=h5py.File(outfile_cc_2,"w")
h5_outcc_3=h5py.File(outfile_cc_3,"w")
h5_outcc_1.create_dataset("angle",data=output_cc_1)
h5_outcc_2.create_dataset("angle",data=output_cc_2)
h5_outcc_3.create_dataset("angle",data=output_cc_3)
h5_outcc_1.close()
h5_outcc_2.close()
h5_outcc_3.close()

print("    No. of volume 1 : "+str(len(output_cc_1))+"\n    No. of volume 2 : "+str(len(output_cc_2))+"\n    No. of volume 3 : "+str(len(output_cc_3))+"\nCompare CC finish")

