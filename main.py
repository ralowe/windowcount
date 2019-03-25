from sys import argv

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def printHelp():
    print("Usage: windowcount -g <chromosome sizes>")
    print("Options:")
    print("\t-g\t\tSet the genome chromosome sizes")
    exit()

def parseChromSizes(filename):
    out=[]
    for k in open(filename):
        c=k.split("\t")[0]
        l=int(k.split("\t")[1])
        out.append([c,l])
    return(out)
    
def main():
    myargs = getopts(argv)
    if '-g' in myargs: 
        ref=myargs["-g"]
    else:
        printHelp()
    cs=parseChromSizes(ref)
        
if __name__=="__main__":
    main()