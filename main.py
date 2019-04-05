from sys import argv
from sys import stderr
import os
import pysam
import time

class BamFile:
    def __init__(self,filename):
        self.itr = pysam.AlignmentFile(filename, "rb").fetch()
        self.readnext()
    
    def currentRead(self):
        return(self.r)
        
    def readnext(self):
        while True:
            self.r=next(self.itr)
            if self.r.next_reference_start>self.r.reference_start:
                self.frag=[self.r.reference_name,self.r.reference_start,self.r.reference_start+self.r.template_length,self.r.query_name]
                break
        
    def getfragment(self):
        return self.frag
        

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
    cs=[]
    ls=[]
    for k in open(filename):
        c=k.split("\t")[0]
        l=int(k.split("\t")[1])
        cs.append(c)
        ls.append(l)
    return(cs,ls)
    
def extractreadcounts(chrom,spos,epos,sample,buffer):
    #Check we havent gone to the next chromosome
    while (sample.frag[1]<epos and chrom==sample.frag[0]):
        buffer.append(sample.frag)
        sample.readnext()
    
    for k in range(len(buffer)-1,-1,-1):
        if buffer[k][2]<spos:
            del buffer[k]
    
    count=0
    if buffer:
        count+=len(buffer)
    return buffer,count
    
def windowcount(chrom,l,size,sep,files,outfile):
    p=0
    readbuffers=[[]]*len(files)
    #Only allow complete windows starting from 0
    while (p+size)<l:
        counts=[]
        for z in range(len(files)):
            readbuffers[z],c1 = extractreadcounts(chrom,p,p+size,files[z],readbuffers[z])
            counts.append(c1)
        if sum(counts)>0:
            string=f"{chrom}\t{p}\t{p+size}"
            for c in counts:
                string+=f"\t{c}"
            outfile.write(f"{string}\n")
        p+=sep
    
def checkchrom(cs,current,files):
    #Move to the right chromosome
    for z in files:
        while(cs.index(z.frag[0])<current):
            z.readnext()

    
def getfiles(direc,outfile):
    files=[]
    header="chrom\tspos\tepos"
    for i in os.listdir(direc):
        if i.find(".bam")!=-1 and i.find(".bai")==-1:
            files.append(BamFile(f"{direc}/{i}"))
            header+=f"\t{i.split('.bam')[0]}"
    outfile.write(f"{header}\n")
    return(files)
        
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
    cs,ls=parseChromSizes(ref)
    outfile=open("test.txt","w")
    files=getfiles(direc,outfile)
    for k in range(len(cs)):
        print("*******************")
        print(k)
        t=time.time()
        checkchrom(cs,k,files)
        windowcount(cs[k],ls[k],w,s,files,outfile)
        stderr.write(f"Done {cs[k]} in {time.time()-t}s\n")
if __name__=="__main__":
    main()