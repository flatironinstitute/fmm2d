
fname="jinjaroot.yaml"

file = open(fname,"w")

file.write("mbhDirectRouts:\n")

outs=["p","g","h"]

for out in outs:
    for i in range(16):
        i1 = i % 2
        i2 = (i // 2) % 2
        i3 = (i // 4) % 2
        i4 = (i // 8) % 2
        ker = ""
        if (i1 == 1): ker += "c"
        if (i2 == 1): ker += "d"
        if (i3 == 1): ker += "q"
        if (i4 == 1): ker += "o"

        if ker != "":
            name = "mbh2d_direct"+ker+out+"_vec"
            file.write("  -\n")
            file.write("    name: " + name + "\n")
            file.write("    ker: " + ker + "\n")
            file.write("    out: " + out + "\n")                        
            
        
