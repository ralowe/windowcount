from sys import argv
import os
import pysam

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def printHelp():
    print("Usage: windowcount -g <chromosome sizes> -d <directory of BAMs>")
    print("Options:")
    print("\t-g\t\tSet the genome chromosome sizes")
    print("\t-d\t\tDirectory to search for BAM files")
    print("\t-w\t\tSize of window [Int: Default 100]")
    print("\t-s\t\tAmount of shift between each window [Int: Default 50]")
    exit()

def parseChromSizes(filename):
    out=[]
    for k in open(filename):
        c=k.split("\t")[0]
        l=int(k.split("\t")[1])
        out.append([c,l])
    return(out)
    
def extractreadcounts(chrom,spos,epos,sample):
    count=0
    for read in sample.fetch(chrom, spos, epos):
        count+=1
    return count
    
def windowcount(chromsizes,size,sep,files):
    l=chromsizes[1]
    p=0
    #Only allow complete windows starting from 0
    while (p+size)<l:
        counts=[]
        for f in files:
            counts.append(extractreadcounts(chromsizes[0],p,p+size,f))
        if sum(counts)>0:
            string=f"{chromsizes[0]}\t{p}\t{p+size}"
            for c in counts:
                string+=f"\t{c}"
            print(string)
        p+=sep
    
def getfiles(direc):
    out=[]
    header="chrom\tspos\tepos"
    for i in os.listdir(direc):
        if i.find(".bam")!=-1 and i.find(".bai")==-1:
            out.append(pysam.AlignmentFile(f"{direc}/{i}", "rb"))
            header+=f"\t{i.split('.bam')[0]}"
    print(header)
    return(out)
        
def main():
    myargs = getopts(argv)
    w=100
    s=50
    if '-g' in myargs: 
        ref=myargs["-g"]
    else:
        printHelp()
    if '-d' in myargs: 
        direc=myargs["-d"]
    else:
        printHelp()
    if '-w' in myargs: 
        w=int(myargs["-w"])
    if '-s' in myargs: 
        s=int(myargs["-s"])
    cs=parseChromSizes(ref)
    files=getfiles(direc)
    for k in cs:
        windowcount(k,w,s,files)
        
if __name__=="__main__":
    main()