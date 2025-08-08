a = "KEVKAAELAAAKEAAKAELKALNLSEGQKDFYIKKINDAKTVEGVKALLEEALKLNDAKKEI"
b = "KEVKEEELLKAKKEAIEELRKLNLSREQKLFYIKKILDAKTVEGVKALLEEALKLNDAKKEI"
c = 0
where_muts = []
for i in range(len(a)):
    if a[i] != b[i]:
        c += 1
        where_muts.append(i + 1)
print(c)
print(where_muts)

wt_marked = "KEVKAAELAAAKXXAXXELXXLNLSXXQXXFYIXKIXDAKTVEGVKALLEEALKLNDAKKEI"

ordered_sites = []
for i, site in enumerate(wt_marked):
    if site == "X":
        ordered_sites.append(i + 1)


print(ordered_sites)