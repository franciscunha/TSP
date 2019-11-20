fin = open("bm.txt","r")
fout = open("results.txt", "w+")

i = 0
flag = 0 # 0 is cost 1 is time
for line in fin.readlines():
    no_slash_n = line.split("\n")
    my_list = no_slash_n[0].split()
    
    if len(my_list) == 1:
        
        if i > 0:
            avg_cost = avg_cost/10
            avg_time = avg_time/10

            fout.write("\nAvg cost: ")
            fout.write(str(avg_cost))
            fout.write("\nAvg time: ")
            fout.write(str(avg_time))
            fout.write("\n\n")

        fout.write(my_list[0])

        avg_cost = 0
        avg_time = 0

    else:
        num = my_list[1]
        
        if flag == 0:
            avg_cost += int(num)
            flag = 1
        elif flag == 1:
            avg_time += float(num)
            flag = 0
        
    i += 1


fin.close()
fout.close()

