import os
for filename in os.listdir(os.getcwd()):
   if filename[-2:] == 'py':
     with open(os.path.join(os.getcwd(), filename), 'r') as f:
      with open(os.path.join(os.getcwd(), filename+'2'), 'w') as f2:
       for line in f:

        if("print(" in line):)
          f2.write(line.replace("print(", "print(").rstrip()+')\n'))
        else:
          f2.write(line)
     os.system('mv '+os.path.join(os.getcwd(), filename)+'2 '+os.path.join(os.getcwd(), filename))
