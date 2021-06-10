from ctypes import c_int, cdll, c_char_p

# go build -buildmode=c-shared -o s1.so main.go reads.go
lib = cdll.LoadLibrary("./s1.so")
lib.Fetch.argtypes = [c_char_p, c_char_p, c_int, c_int]
lib.Fetch.restype = c_char_p
result = lib.Fetch(
    "/mnt/raid64/ATS/Personal/zhangyiming/bams/NCovM1.bam".encode("utf-8"), 
    "1".encode("utf-8"), 
    1000, 
    1200
)


print(result.decode("utf-8"))