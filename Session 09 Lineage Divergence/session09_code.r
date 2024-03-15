gl.set.verbosity(3)

gl <- readRDS("./data/Example_12.1_SNP.Rdata")

gl

nLoc(gl)

nInd(gl)

nPop(gl)

indNames(gl)[1:10]

popNames(gl)

table(pop(gl))

gl@other$loc.metrics$TrimmedSequence[1:10]

gl.report.taglength(gl)

gl@other$loc.metrics$SnpPosition[1:10]
                              
gl@position[1:10]

gl@other$loc.metrics$SNP[1:10]

gl.report.bases(gl)

gl2fasta(gl,outfile="Example1.fas",outpath = './data/')

gl <- gl.filter.taglength(gl, lower=40)

d <- gl.dist.phylo(gl,subst.model="F81")
gl.tree.fitch(D=d,x=gl,bstrap=100, )


gl.tree.nj(d)


gl <- gl.load("Example_SNP.Rdata")

# Consider if you want to filter on taglength, say
gl <- gl.filter.taglength(gl,lower=69, verbose=3)


dd <-gl.dist.phylo(gl,by.pop=TRUE, subst.model="F81",verbose=3)

as.matrix(dd)

gl.tree.fitch(D=dd,x=gl, phylip.path="D:/workspace/R/phylip-3.695/exe",
              outgroup="macquarii", verbose=3)


gl.tree.fitch(D=dd,x=gl, phylip.path="D:/workspace/R/phylip-3.695/exe",
              outgroup="macquarii", bstrap=100, verbose=3)

