
#import numpy as np
path = './sequence.fasta'
f = open(path)
s1 = f.read()
data = "".join(s1.split("\n")[1:])
#sequence = 'CGCTAG'
x = y = 0.0
x_array = [0]
y_array = [0]


for i in data:
    if i == 'A':
        x = (-1 + x)/2
        y = (1 + y)/2
        x_array.append(x)
        y_array.append(y)
    elif i == 'G':
        x = (1 + x)/2
        y = (-1 + y)/2
        x_array.append(x)
        y_array.append(y)
    elif i == 'C':
        x = (-1 + x)/2
        y = (-1 + y)/2
        x_array.append(x)
        y_array.append(y)
    elif i == 'T':
        x = (1 + x)/2
        y = (1 + y)/2
        x_array.append(x)
        y_array.append(y)
    else:
        pass

#print x_array[-2]
x = y = 0.0
#y1 = y_array[-1] * 2
#x1 = y_array[-1] * 2
x_a = []
y_a = []
pre_sequence = ''
for i in range(len(x_array)):
    x = x_array[-i-1] * 2
    x_a.append(x)

for i in range(len(y_array)):
    y = y_array[-i-1] * 2
    y_a.append(y)


pre_x = x_a[::-1]
pre_y = y_a[::-1]
for j,k in zip(pre_x,pre_y):
    if j < 0 and k > 0:
        pre_sequence += 'A'
    elif j > 0 and k > 0:
        pre_sequence += 'T'
    elif j < 0 and k < 0:
        pre_sequence += 'C'
    elif j > 0 and k < 0:
        pre_sequence += 'G'
    else:
        pass


#print 'Data value: ' + str(data[53668:53672])
#print 'y_array vale: ' + str(y_array[53668:53672])

if pre_sequence == data:
    print 'Success'
count = 0
wrong_index = []
for i in range( 0, len(pre_sequence)):
    count += 1
    if pre_sequence[i] != data[i]:
        wrong_index.append(count)
print wrong_index 
