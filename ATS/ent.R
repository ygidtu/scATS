fs = list.files("testfile/", pattern = "*", full.name = T)

for (f in fs) {
    Z = read.table(f)
    LZ = Z
    LZ[Z!=0] = log(Z[Z!=0])
    entropy = -1 * Z * LZ
    sum(entropy)
}
