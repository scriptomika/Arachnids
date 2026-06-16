#run Dollo on server
# input: Rdata with
# -character table (species rows, characters columns), 
# -rooted tree (class phylo)
# -vector of nodenames
CORES<-23
# (plotting commands suppressed)
require(phytools)
library(phangorn)

load("og4dollo.Rdata") #loads ogdf, tr1, tr2 nodenames1, nodenames2
#TESTING tree hypothesis 1
tree<-tr1
nodenames<-nodenames1

runDollo<-function(og,sigOGs,tree){
  # og<-"OG0000963" #use for testing 
  char<-sigOGs[tree$tip.label,og]
  
  #plot(tree); axisPhylo()
  #title(colnames(sigOGs)[og])
  #tiplabels(tip=which(tree$tip.label %in% tree$tip.label[char==0]),col=absent.col,pch=16,cex=1,frame="none")
  #tiplabels(tip=which(tree$tip.label %in% tree$tip.label[char>0]),col=present.col,pch=16,cex=1,frame="none")
  
  #Dollo step1: establish list of tips with Present state
  present.tips<-which(tree$tip.label %in% names(which(char>0))) #char is already matching tree tip.label order
  anc.nodes<-c()
  
  
  for (dec in 1:length(present.tips)){
    #Dollo step 2: determine youngest nodes subtending Present tips
    parentnode<-tree$edge[which(tree$edge[,2]==present.tips[dec]),1]
    anc.nodes<-c(anc.nodes,parentnode)
  }
  anc.nodes<-unique(anc.nodes)
  #nodelabels(node=anc.nodes,pch=16,col=present.col)
  
  # Dollo step 3: For the remaining, yet-undetermined nodes: save as Present IF there are both older AND younger nodes saved as Present. Any MRCA between Present nodes must also be a Present.
  
  # capture descendents (nodes+tips) for each node
  descs<-phangorn::allDescendants(x=tree)
  
  # create list of undetermined nodes by using the parent nodes of just-determined (step2) Present nodes to subset the whole tree's internal nodes and retain only those internal nodes representing a MRCA to Present nodes
  intnodes<-unique(tree$edge[,1])
  nextnodeup<-unique(tree$edge[which(tree$edge[,2] %in% anc.nodes),1])
  rem_intnodes<-intnodes[intnodes>=min(nextnodeup)]
  
  # ensure root is excluded for now (later, explicitly check if outgroup descendents are Present before saving root as Present)
  root.node<-Ntip(tree)+1
  rem_intnodes<-rem_intnodes[!(rem_intnodes %in% root.node)]
  
  #for the node preceeding the last certain Present node, check if it has one absent descendent tip, if so exclude this node from next step
  if(all(Descendants(tree,min(rem_intnodes),"children") %in% c(present.tips,anc.nodes))==FALSE)
  {rem_intnodes<-rem_intnodes[!(rem_intnodes %in% min(rem_intnodes))]
  }
  
  
  # save these undetermined nodes as Present if any of their descendents are Present
  for (n in unique(rem_intnodes)){
    if (any(present.tips %in% unlist(descs[n]))){ anc.nodes<-c(anc.nodes,n)}
  }
  anc.nodes<-unique(anc.nodes)
  #nodelabels(node=anc.nodes,pch=16,col=present.col)
  
  # Dollo step 4: check if the deepest MRCAs between Present descendents hasn't been captured as Present
  ti<-c(unlist(descs[min(anc.nodes)]),unlist(descs[max(anc.nodes)])) 
  tips<-tree$tip.label[ti [ti<= Ntip(tree)]]
  mrca1<-getMRCA(tree,tips)
  if (any(present.tips %in% unlist(descs[mrca1]))){ anc.nodes<-c(anc.nodes,mrca1)}
  #nodelabels(node=anc.nodes,pch=16,col=present.col)
  
  #if there are intermediate nodes between mrca-present and presnet-nodes, add to present list
  #overlap of getDescendents and getAncesotrs(of present nodes)
  for (n in anc.nodes){
    ancs.tmp<-phangorn::Ancestors(x=tree,node=n)
    anc.nodes<-c(anc.nodes,ancs.tmp [ancs.tmp %in% unlist(descs[mrca1])])
  }
  anc.nodes<-unique(anc.nodes)
  #nodelabels(node=anc.nodes,pch=16,col=present.col)
  
  
  # Dollo step 5: if both of the desendent nodes of the root are in the Present node list, add root to list
  rootchildren<-tree$edge[tree$edge[,1]==root.node,2]
  if (any(present.tips %in% unlist(descs[rootchildren[21]])) & any(present.tips %in% unlist(descs[rootchildren[2]]))){ anc.nodes<-c(anc.nodes,root.node)}
  
  # save results in table
  presnodes<-anc.nodes-Ntip(tree)
  absnodes<-unique(sort(tree$edge[,1]))[!(unique(sort(tree$edge[,1])) %in% anc.nodes)]-Ntip(tree)
 # anc.states[ presnodes,og]<-"present"
 # anc.states[absnodes,og]<-"absent"
  
  #instead create one vector Nnode long combining pres and abs
  pn<-rep(1,length(presnodes)); pn<-cbind(pn,presnodes)
  an<-rep(0,length(absnodes)); an<-cbind(an,absnodes)
  nodestates<-rbind(pn,an);nstates<-nodestates[,"pn"][order(nodestates[,"presnodes"])]
return(nstates)
  
  }

#ogdf_test<-ogdf[,1:10]
#listnodestates<-lapply(colnames(ogdf_test),FUN=function (x) runDollo(x,ogdf_test,tree) )
#names(listnodestates)<-colnames(ogdf_test)

if (!file.exists("dolloASR_results_H1.Rds")){
cat("Starting Dollo ASR predictions...\n")
  st<-Sys.time()
library(parallel)                  
listnodestates<-mclapply(colnames(ogdf),function (x) runDollo(x,ogdf,tree), mc.cores=CORES)
names(listnodestates)<-colnames(ogdf)
en<-Sys.time()
dollo.anc.states<-data.frame(listnodestates)
rownames(dollo.anc.states)<-nodenames

#dollo.anc.states[dollo.anc.states=="present"]<-1
#dollo.anc.states[dollo.anc.states=="absent"]<-0
#dollo.anc.states<-as.data.frame(lapply(dollo.anc.states, as.numeric))

#save object
save(dollo.anc.states,file="dolloASR_results_H1.Rds")
cat("Finished Dollo ASR results table (file 1/3)")
print(en-st)
}else{load("dolloASR_results_H1.Rds")}

#tabulate shifts
if (!file.exists("loss_gain_pres_abs_H1.Rds")){
root.node<-Ntip(tree)+1

rownames(dollo.anc.states)<-seq(from=(Ntip(tree)+1),to=(Nnode(tree)+Ntip(tree)))
char<-ogdf[tree$tip.label,]
alledges<-rbind(char,dollo.anc.states)
rownames(alledges)<-seq(from=1,to=nrow(tree$edge)+1)

assess_shifts<-function(edgerow,alledges,og){ 
  ancnode<-edgerow[1];# print(ancnode)
  decnode<-edgerow[2];#print(decnode)
  if (alledges[ancnode,og]>0) {
     r<-ifelse(alledges[decnode,og]==0,"l","p")}
  else{
    r<-ifelse(alledges[decnode,og]==0,"a","g")}
  return(c(decnode,r))}


#lossgain.table<-matrix(NA,nrow=nrow(tree$edge),ncol=ncol(dollo.anc.states),dimnames = list(rownames(tree$edge),colnames(dollo.anc.states)))
#for (og in colnames(alledges)){
#  outall<-apply(tree$edge,1,FUN=function(x) assess_shifts(x,alledges,og))
#  lossgain.table[,og]<-outall[2,]
#}
cat("Starting tabulation of gains and losses between nodes...\n")
st<-Sys.time()
library(parallel)
outlist<-mclapply(colnames(alledges), FUN= function(i) apply(tree$edge,1,FUN=function(x) assess_shifts(x,alledges,i)),mc.cores=CORES)
names(outlist)<-colnames(alledges)

node.order<-order(as.numeric(unlist(outlist[[1]])[1,]))
outlist2<-lapply(outlist,FUN= function(x) unlist(x)[2,][node.order])
lossgain.table<-data.frame(outlist2)

en<-Sys.time()
# root (115 has no corresponding rows)
rootstates<-dollo.anc.states[1,] %>% slice(1) %>% unlist() %>% dplyr::recode(.,`0`="a",`1`="p")
lossgain.table<-rbind.data.frame(lossgain.table[1:Ntip(tree),],rootstates,lossgain.table[(Ntip(tree)+1):nrow(lossgain.table),])
rownames(lossgain.table)<-1:nrow(lossgain.table)
#now root is in table, but no gain/loss shown as would require comparison with earlier node


#extract nodes for plotting. action on the decendent node, since it is the consequnece of gain/loss on subtending branch.
#gain loss table1 shows outcome in dec node for each row of anc-dec pairs in tree$edge structure, same order
# e.g. 
# tree edge    table 227 rows == number of nodes and tip (all descendents)
# 115 116       a
# 116 1         a
# 116 2         a
# 115 117       g
# 
# root row will only show p or a as no gain/loss may be infered
print("Finished loss/gain table")
print(en-st)
save(lossgain.table,file="loss_gain_pres_abs_H1.Rds")
cat("Finished gain/loss table (file 2/3)")
}else{load("loss_gain_pres_abs_H1.Rds")}

cat("Preparing final summary table...\n")
  gainloss.dollo.mat<-t(apply(lossgain.table,MARGIN=1, function(x) table(factor(x, levels=c("a","g","l","p")))))
  
  gainloss.dollo.mat.nodes<-gainloss.dollo.mat[seq(from=(Ntip(tree)+1),to=Nnode(tree)+Ntip(tree)-1),]
  
  #root is missing in table, since gain/loss require comparison with previous node
  rownames(gainloss.dollo.mat.nodes)<-nodenames[-1] #drop root
  save(gainloss.dollo.mat.nodes,file="gainloss.summary_H1.Rds")

  #save lists of OGs for each node
  presOGlist<-apply(dollo.anc.states,1,FUN=function(x) names(x[which(x==1)]))
  absOGlist<-apply(dollo.anc.states,1,FUN=function(x) names(x[which(x==0)]))
  
  lossgain.t2<-apply(lossgain.table,2,as.character)
  #label rows with nodenames and tip names
  rownames(lossgain.t2)<-c(tree$tip.label,as.character(nodenames[-1]))
  
  gainOGlist<-apply(lossgain.t2,1,FUN=function(x) names(x[which(x=="g")]))
  presOGlist<-apply(lossgain.t2,1,FUN=function(x) names(x[which(x=="p")]))
  lostOGlist<-apply(lossgain.t2,1,FUN=function(x) names(x[which(x=="l")]))
  absentOGlist<-apply(lossgain.t2,1,FUN=function(x) names(x[which(x=="a")]))
  
  #list of genes present at root:
  root_present<-presOGlist[1]
  
  save(root_present,gainOGlist,lostOGlist,presOGlist,absentOGlist,file="OGlists_by_node_tip_H1.Rds")
cat("Saved lists of OGs by node (file 3/3)")
