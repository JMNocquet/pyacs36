

def f():
    global s
    s = 10.
    
    return
    
def g():
    global s
    f()
    print(s)
    s = 13.

def h():
    i = 12.
    return

g()
print(s)
