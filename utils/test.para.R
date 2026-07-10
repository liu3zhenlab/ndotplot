datafile="test/out/02_delta.txt"
lend.turnoff=F
refname="ref"
qryname="asm"
xlabel.rm = "chr"
ylabel.rm = "tig0+"
line.width.factor=3.5
#outpdf

rorder <- c("ref")
qorder <- c("qry")


tstr <- "a, b,c"
tvec <- strsplit(tstr, "[,;]")[[1]]
tvec
length(tvec)
