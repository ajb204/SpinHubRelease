
import os
filey =os.listdir('src')
for file in filey:
    inny=open('srcCopy/'+file)
    outy=open('tmp','w')
    for line in inny.readlines():
        line=line.replace('prodNew(','propagate(')
        line=line.replace('Sigma,Sigma);','Sigma);')
        line=line.replace('Sigma, Sigma);','Sigma);')
        outy.write(line)
    outy.close()
    os.system('mv tmp src/'+file)
    inny.close()
