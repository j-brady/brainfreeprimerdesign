#!/Users/jacobbrady/virtual_envs/py35/bin/python
import sys
print(sys.argv[1:])
seq = "".join(sys.argv[1:])
seq_str, melt_temp, melt_temp_list = "", 0, ""
for i in seq.upper().strip():
    if (i == "A") or (i == "T"):
        seq_str += "%4s"%i
        melt_temp += 2
        melt_temp_list += "%4d"%melt_temp
    elif (i == "G") or (i == "C"):
        seq_str += "%4s"%i
        melt_temp += 4
        melt_temp_list += "%4d"%melt_temp
print(seq_str)
print(melt_temp_list)
