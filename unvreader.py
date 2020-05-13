f = open('Mesh.unv','r')
lines = f.readlines()
lines = lines.split()
lines = lines[22:]


#read 2411 nodes
for i,l in enumerate(lines):
    if len(l)==1 and float(l[0])==-1.:
        break
    