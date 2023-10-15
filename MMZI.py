from inspect import signature, isfunction
import sympy as sp
import numpy as np
import math


def matrica_prijelaza(baza1,baza2,klasa='sympy'):
    """Daje matricu prijelaza iz baze baza1 u bazu baza2.
    Opcija klasa moze imati dvije vrijednosti 'sympy' ili 'numpy' ovisno
    koju klasu matrice zelimo dobiti na izlazu.
    
    --- PRIMJER ---
     
     B1=[(1,0,0),(1,1,0),(3,4,5)]
     B2=[(3,1,-1),(1,-1,-1),(0,0,2)]
     matrica_prijelaza(B1,B2)
     matrica_prijelaza(B1,B2,'numpy')
    """
    if klasa == 'sympy':
        m1 = sp.Matrix(baza1).T
        if m1.det() == 0: return "Error: prvi skup nije baza"
        m2 = sp.Matrix(baza2).T
        if m2.det() == 0: return "Error: drugi skup nije baza"
        if m1.shape != m2.shape: return "Error: baze nisu iz istog vektorskog prostora"
        return m1.inv()*m2
    elif klasa == 'numpy':
        m1 = np.transpose(np.array(baza1))
        if np.abs(np.linalg.det(m1)) < 1e-15: return "Error: prvi skup nije baza"
        m2 = np.transpose(np.array(baza2))
        if np.abs(np.linalg.det(m2)) < 1e-15: return "Error: drugi skup nije baza"
        if m1.shape != m2.shape: return "Error: baze nisu iz istog vektorskog prostora"
        return np.matmul(np.linalg.inv(m1), m2)


def operator_kan(f,klasa='sympy'):
    """Daje matricni zapis linearnog operatora f u paru kanonskih baza.
    Linearni operator mora biti prethodno definiran kao funkcija ili moze biti
    definiran na ulazu kao lambda funkcija.
    
    Opcija klasa moze imati dvije vrijednosti 'sympy' ili 'numpy' ovisno
    koju klasu matrice zelimo dobiti na izlazu.
    
    --- PRIMJER ---
     
     def f(x,y):
         return (x+y,x-2*y,3*y)
     
     operator_kan(f)
     
     operator_kan(lambda x,y:(x+y,x-2*y,3*y))
    """     
    n = len(signature(f).parameters)
    if klasa == 'sympy':
        identiteta = sp.eye(n)
        kan = [list(identiteta.col(i)) for i in range(identiteta.cols)]
        slika = list(map(lambda t: f(*t), kan))
        matrica = sp.Matrix(slika).T
    elif klasa == 'numpy':
        kan = np.eye(n)
        slika = np.apply_along_axis(lambda t: f(*t),1,kan)
        matrica = np.transpose(slika)
    return matrica


def operator_kanbaza(f,baza1,baza2=None,klasa='sympy'):
    """Daje matricni zapis linearnog operatora f u paru baza (baza1,baza2).
    Ako se baza2 eksplicitno ne navede, podrazumijeva se da je baza2=baza1.
    
    Linearni operator f moze biti:
      * prethodno definiran kao funkcija 
      * definiran na ulazu kao lambda funkcija
      * definiran svojim matricnim zapisom u paru kanonskih baza, a pritom
        matrica mora pripadati klasi koja je specificirana u varijabli 'klasa'

    Opcija klasa moze imati dvije vrijednosti 'sympy' ili 'numpy' ovisno
    koju klasu matrice zelimo dobiti na izlazu.
    
    --- PRIMJER ---
     
     def f(x,y):
         return (x+y,x-2*y,3*y)
     
     operator_kanbaza(f,[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)])

     operator_kanbaza(lambda x,y:(x+y,x-2*y,3*y),[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)])

     operator_kanbaza(sp.Matrix([[1,1],[1,-2],[0,3]]),[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)],klasa='sympy')
     operator_kanbaza(np.array([[1,1],[1,-2],[0,3]]),[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)],klasa='numpy')
    """ 
    if isfunction(f):
        #operator f je zadan formulom
        M = operator_kan(f,klasa)
    else:
        #operator f je zadan matricom u paru kanonskih baza
        M = f
    if klasa == 'sympy':
        S = sp.Matrix(baza1).T
        if S.det() == 0: return "Error: prvi skup nije baza"
        if baza2 == None:
            T = sp.Matrix(baza1).T
        else:
            T = sp.Matrix(baza2).T
            if T.det() == 0: return "Error: drugi skup nije baza"
        return T.inv()*M*S
    elif klasa == 'numpy':
        S = np.transpose(np.array(baza1))
        if np.abs(np.linalg.det(S)) < 1e-15: return "Error: prvi skup nije baza"
        if baza2 == None:
            T = np.transpose(np.array(baza1))
        else:
            T = np.transpose(np.array(baza2))
            if np.abs(np.linalg.det(T)) < 1e-15: return "Error: drugi skup nije baza"
        return np.matmul( np.matmul(np.linalg.inv(T), M), S)


def operator_baza(f,bazaD1,bazaD2,bazaK1=None,bazaK2=None,klasa='sympy'):
    """ Ova funkcija je poopcenje funkcije operator_kanbaza. Ona daje
    matricni prikaz linearnog operatora f u paru baza (bazaD2,bazaK2)
    ako je on zadan u paru baza (bazaD1,bazaK1).
    Ako bazaK1 nije eksplicitno navedena, podrazumijeva se da je bazaK1=bazaD1.
    Ako bazaK2 nije eksplicitno navedena, podrazumijeva se da je bazaK2=bazaD2.
    
    Linearni operator f moze biti:
      * prethodno definiran kao funkcija u paru baza (bazaD1,bazaK1)
      * definiran na ulazu kao lambda funkcija u paru baza (bazaD1,bazaK1)
      * definiran svojim matricnim zapisom u paru baza (bazaD1,bazaK1), 
        pri cemu matrica mora pripadati klasi navedenoj u varijabli 'klasa'
    
    --- 1. PRIMJER ---
    
    Neka je [[2,0],[-1,0]] matrica operatora u paru baza {(1,1),(1,0)} i {(1,0),(0,1)}.
    Trazimo matricu tog operatora u bazi B={(1,0),(1,1)}.
    Dakle, u ovom slucaju je bazaD1={(1,1),(1,0)}, bazaD2=B, bazaK1={(1,0),(0,1)}, bazaK2=B.
    U ovom slucaju je bazaK2=bazaD2 pa ju nema potrebe dolje pisati na zadnjem mjestu.
    
    operator_baza(sp.Matrix([[2,0],[-1,0]]), [(1,1),(1,0)], [(1,0),(1,1)], [(1,0),(0,1)])
    operator_baza(np.array([[2,0],[-1,0]]), bazaD1=[(1,1),(1,0)], bazaD2=[(1,0),(1,1)], bazaK1=[(1,0),(0,1)], klasa='numpy')
    
    --- 2. PRIMJER ---
    
    Neka je [[1,0,2],[1,1,1],[1,0,1]] matrica linearnog operatora u bazi {(1,1,1),(1,2,1),(0,1,2)}.
    Trazimo matricu tog operatora u bazi {(1,0,1),(-1,0,1),(0,2,0)}.
    U ovom slucaju je bazaD1={(1,1,1),(1,2,1),(0,1,2)}, bazaD2={(1,0,1),(-1,0,1),(0,2,0)}, bazaK1=bazaD1, bazaK2=bazaD2.
    Kako je bazaK1=bazaD1, bazaK2=bazaD2, nema potrebe posebno pisati baze bazaK1, bazaK2 kod poziva funkcije.
    
    operator_baza(sp.Matrix([[1,0,2],[1,1,1],[1,0,1]]), [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)])
    operator_baza(np.array([[1,0,2],[1,1,1],[1,0,1]]), [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)], klasa='numpy')
    
    Nije greska ako ih i posebno naglasimo.
    
    operator_baza(sp.Matrix([[1,0,2],[1,1,1],[1,0,1]]), [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)], [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)])
    
    ili da nam bude preglednije
    
    operator_baza(sp.Matrix([[1,0,2],[1,1,1],[1,0,1]]), bazaD1=[(1,1,1),(1,2,1),(0,1,2)], bazaD2=[(1,0,1),(-1,0,1),(0,2,0)])
    operator_baza(sp.Matrix([[1,0,2],[1,1,1],[1,0,1]]), bazaD1=[(1,1,1),(1,2,1),(0,1,2)], bazaD2=[(1,0,1),(-1,0,1),(0,2,0)], 
                  bazaK1=[(1,1,1),(1,2,1),(0,1,2)], bazaK2=[(1,0,1),(-1,0,1),(0,2,0)])
    operator_baza(np.array([[1,0,2],[1,1,1],[1,0,1]]), bazaD1=[(1,1,1),(1,2,1),(0,1,2)], bazaD2=[(1,0,1),(-1,0,1),(0,2,0)], 
                  bazaK1=[(1,1,1),(1,2,1),(0,1,2)], bazaK2=[(1,0,1),(-1,0,1),(0,2,0)], klasa='numpy')
    """     
    if isfunction(f):
        #operator f je zadan formulom
        M = operator_kan(f,klasa)
    else:
        #operator f je zadan matricom u paru baza (bazaD1,bazaK1)
        M = f
    if bazaK1 == None: bazaK1 = bazaD1
    if bazaK2 == None: bazaK2 = bazaD2
    S = matrica_prijelaza(bazaD1,bazaD2,klasa)
    T = matrica_prijelaza(bazaK2,bazaK1,klasa)
    if klasa == 'sympy':
        return T*M*S
    elif klasa == 'numpy':
        return np.matmul(np.matmul(T,M),S)


def razvoj_der(pol,br):
    """ Funkcija daje koeficijente u razvoju polinoma pol po potencijama od x-br.
    Varijabla pol mora biti klase poly1d. Na izlazu se daje lista koeficijenata
    pri cemu je element u pripadnoj listi s indeksom k jednak koeficijentu koji 
    stoji uz potenciju (x-br)^k. Implementacija koristi pristup preko derivacija.

    --- PRIMJER ---(numpy modul je ucitan kao: import numpy as np)

    Razvoj polinoma p(x)=-3x^4+5x^2+x-8 po potencijama od x+5.
    p=np.poly1d([-4,0,5,1,-8])
    razvoj_der(p,-5)
    """
    razvoj=[pol.deriv(k)(br)/math.factorial(k) for k in range(pol.o+1)]
    return razvoj


def razvoj_horner(pol,br):
    """ Funkcija daje koeficijente u razvoju polinoma pol po potencijama od x-br.
    Varijabla pol mora biti klase poly1d. Na izlazu se daje lista koeficijenata
    pri cemu je element u pripadnoj listi s indeksom k jednak koeficijentu koji 
    stoji uz potenciju (x-br)^k. Implementacija koristi Hornerov algoritam.

    --- PRIMJER ---(numpy modul je ucitan kao: import numpy as np)

    Razvoj polinoma p(x)=-3x^4+5x^2+x-8 po potencijama od x+5.
    p=np.poly1d([-4,0,5,1,-8])
    razvoj_horner(p,-5)
    """
    kv=pol
    lin=np.poly1d([1,-br])
    raz=[]
    for k in range(pol.o+1):
        kv,ost=kv/lin
        raz.append(ost[0])
    return raz


def Ferrari_rezolventa(koef):
    """Daje Ferrarijevu rezolventu algebarske jednadzbe 4. reda. Varijabla koef
    mora biti lista od 5 brojeva koji predstavljaju koeficijente jednadžbe od
    najvece prema najmanjoj potenciji. Ferrarijeva rezolventa je na izlazu 
    reprezentirana preko numpy klase poly1d.

    ---PRIMJER---
    Ferrari_rezolventa([-3,4,0,6.78,-2])
    """

    if len(koef) != 5:
        return "Error: mora biti 5 koeficijenata"
    if koef[0] == 0:
        return "Error: prvi koeficijent ne smije biti 0"
    koef = np.array(koef)
    if koef[0] != 1:
        koef = koef / koef[0]
    b1 = 8*koef[4] - 2*koef[1]*koef[3]
    b0 = koef[3]**2 + koef[1]**2*koef[4] - 4*koef[2]*koef[4]
    return np.poly1d([-8,4*koef[2],b1,b0])
