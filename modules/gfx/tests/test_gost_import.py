import struct

with open("test.gost","rb") as gf:

  header = gf.read(6)
  if str(header[0:4]) != "GOST":
    raise RuntimeError("file format mismatch")
  raw = gf.read(16)

  while raw:
    (type, subtype,size) = struct.unpack("iiL",raw)
    print("found type=%d, subtype=%d and blocksize=%u"%(type,subtype,size))
    if size>0:
      data = gf.read(size)
    raw = gf.read(16)
