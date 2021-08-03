load('../filtered_nets.RData')
source('../HirNet_function.R')
data=AD.Ex.filtered[,1:2]
cal_hier_score26(data=data, kmax=10000, ptim=100, anneal.coeff=1e-6,  myoutf = "AD_Ex.result.txt")

data=Ctrl.Ex.filtered[,1:2]
cal_hier_score26(data=data, kmax=10000, ptim=100, anneal.coeff=1e-6,  myoutf = "Ctrl_Ex.result.txt")
