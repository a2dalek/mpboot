
x = [0] * 101
y = [0] * 101
for i in range(1,75):
    f = open("data." + str(i) + ".orig.boot.treefile.rfinfo", "r")
    content = f.read()
    content1 = content.split()
    for num_string in content1:
        num = int(num_string)
        if (num >= 0):
            x[num] += 1
        y[abs(num)] += 1

for i in range(0,101):
    if (y[i] == 0):
        print(str(i) +  " 0 0 0")
    else:
        print(str(i) +  " " + str(x[i]) + " " + str(y[i]) + " " + str(float(x[i])/y[i]))