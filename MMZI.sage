MMZI_naredbe={0:['DMS','bazicna_rjesenja','linear_solve','solve_mat'],
             1:['ProjekcijaVektorVektor','ProjekcijaVektorRavnina','KutVektora','SP','VP','MP'],
             2:['dTocka','Omjer','Ravnina3Tocke','Ravnina2Vektora','RavninaNormala','dTockaRavnina','KutRavnina',
                'RavninaNormalniOblik','Pravac2Tocke','PresjekRavnina','KutPravaca','KutRavninaPravac','dTockaPravac',
                'PresjekPravaca','PresjekPravacRavnina','ZajednickaNormala','UdaljenostPravaca','ProjekcijaTockaRavnina',
                'SimetricnaTockaRavnina','ProjekcijaPravacRavnina','SimetricniPravacRavnina','ProjekcijaTockaPravac',
                'NormalaTockaPravac','SimetricnaTockaPravac','PovrsinaTrokuta','OpsegTrokuta'],
             3:['matrica_prijelaza'],
             4:['operator_kan','operator_kanbaza','operator_baza'],
             5:['rucni_Horner','razvoj_Horner','polinomi_gcd','kompleksni_korijen','Ferrari_rezolventa','Cardan'],
             6:[],
             7:[]}


print("Dodatne funkcije za predmet 'Matematicke metode za informaticare'.")
print("MMZI_naredbe je rjecnik u kojemu je po poglavljima dan popis svih dodatnih naredbi koje trenutno postoje.")


def pojednostavni(x):
    try:
        x=x.simplify_full()
    except:
        pass
    return x

#### 0. poglavlje ####

def DMS(kut,mjera='stupanj'):
    """Daje kut u stupnjevima, minuta, sekundama i stotinkama sekundi.
       Opcija mjera='stupanj' (default) govori da je pocetni kut unesen u stupnjevima.
       Opcija mjera='radijan' govori da je pocetni kut unesen u radijanima.
       
       -- PRIMJERI --
       DMS(40.5678)
       DMS(2.56789,'radijan')"""
    if mjera=='stupanj':
        kut2=kut
    elif mjera=='radijan':
        kut2=(180/pi)*kut
    else:
        return "Error: opcija mjera mora biti 'stupanj' ili 'radijan'"
    stupanj=floor(kut2)
    minuta=(kut2-stupanj)*60
    sekunda=(minuta-floor(minuta))*60
    minuta=floor(minuta)
    sekunda=round(sekunda,2)
    if floor(sekunda)==60:
        sekunda=0
        if minuta<59:
            minuta += 1
        else:
            minuta=0
            stupanj += 1
    rez="{}* {}' {:.2f}''".format(stupanj,minuta,sekunda)
    return rez


def bazicna_rjesenja(matrica,koeficijenti,nepoznanice):
    """Daje sva bazicna rjesenja sustava linearnih jednadzbi.
       Matrica sustava moze biti dvodimenzionalna python lista ili pak moze biti tipa 'matrix'.
       Matrica slobodnih koeficijenata moze biti jednodimenzionalna python lista ili pak moze biti tipa 'matrix'.
       Matrica nepoznanica mora biti jednodimenzionalna python lista.

       --- PRIMJER ---
         
          var('x y z u')
          A=[[1,1,2,3],[1,-4,2,-3],[0,7,2,5],[1,3,4,2]]
          A1=matrix([[1,1,2,3],[1,-4,2,-3],[0,7,2,5],[1,3,4,2]])
          B=[6,-4,12,8]
          B1=matrix(4,1,[6,-4,12,8])
          bazicna_rjesenja(A,B,[x,y,z,u])
          bazicna_rjesenja(A1,B1,[x,y,z,u])
          bazicna_rjesenja(A1,B,[x,y,z,u])
    """
    if type(matrica)==list: matrica=matrix(matrica)
    if type(koeficijenti)!=list: koeficijenti=koeficijenti.list()
    if matrix(matrica.columns()+[tuple(koeficijenti)]).rank()!=matrica.rank():
        return "Error: sustav je kontradiktoran"
    brojParametara=len(nepoznanice)-matrica.rank()
    if brojParametara==0:
        return "Error: sustav ima jedinstveno rjesenje"
    X=matrix(len(nepoznanice),1,nepoznanice)
    lijevo=matrica*X
    desno=matrix(len(koeficijenti),1,koeficijenti)
    jednadzbe=[]
    for i in range(lijevo.nrows()): 
         for j in range(lijevo.ncols()): 
             jednadzbe.append(lijevo[i,j]==desno[i,j])
    parametri=Combinations(zip(nepoznanice,[0]*len(nepoznanice)),brojParametara)
    for t in parametri:
        var1=list(map(lambda q: q[0], t))
        var2=list(filter(lambda q: not(q in var1), nepoznanice))
        rjesenje=solve(list(map(lambda q: q.substitute(dict(t)), jednadzbe)),*var2)
        if rjesenje!=[]:
            rjesenje=list(map(lambda q: q[0]==q[1],t))+rjesenje[0]
            print(rjesenje)
    return None


def linear_solve(lijevo,desno,varijable,parametri=None,rjecnik=False):
    """Rjesava sustav linearnih jednadzbi zapisan u matricnom obliku.
    
    lijevo -> matrica sustava
    desno  -> matrica slobodnih koeficijenata
    varijable -> python lista nepoznanica
    parametri -> redoslijed nepoznanica kojim ce se birati parametri
                 (parametri se biraju s kraja liste prema pocetku).
                 Ako nije eksplicitno navedeno, tada je parametri=varijable.
    rjecnik -> True ili False u ovisnosti zelimo li na izlazu rjesenje
               u obliku python rjecnika ili ne zelimo.
    
    --- PRIMJER ---
     
     x1,x2,x3,x4=var('x1 x2 x3 x4')
     A=matrix(3,4,[2,3,2,6,-2,3,-6,12,2,6,0,15])
     B=matrix(3,1,[1,-19,-8])
     
     linear_solve(A,B,[x1,x2,x3,x4],[x1,x2,x4,x3])
     linear_solve(A,B,[x1,x2,x3,x4])
     linear_solve(A,B,[x1,x2,x3,x4],[x2,x4,x3,x1])
     linear_solve(A,B,[x1,x2,x3,x4],[x2,x4,x3,x1],rjecnik=True)
    """
    jednadzbe=[]
    if parametri==None: parametri=varijable
    lijeva_strana=lijevo*matrix([[t] for t in varijable])
    for i in range(lijeva_strana.nrows()):
        for j in range(lijeva_strana.ncols()):
            jednadzbe.append(lijeva_strana[i,j]==desno[i,j])
    rj=solve(jednadzbe,*parametri,solution_dict=rjecnik)
    return rj


def solve_mat(lijevo,desno,varijable,rjecnik=False):
    """Rjesava matricne jednadzbe tako da se ne moraju raspisivati u sustave linearnih jednadzbi.
    
    lijevo -> lijeva strana matricne jednadzbe
    desno  -> desna strana matricne jednadzbe
    varijable -> nepoznanice u matricnoj jednadzbi (python lista ili tuple)
    rjecnik -> True ili False u ovisnosti zelimo li na izlazu rjesenje
               u obliku python rjecnika ili ne zelimo
    
    --- PRIMJER ---
     
     var('a b c d')
     A=matrix([[1,1],[0,1]])
     X=matrix([[a,b],[c,d]])
     
     solve_mat(A*X,X*A,(a,b,c,d))
     solve_mat(A*X,X*A,(a,b,c,d),rjecnik=True)
     
     X.subs(d=a,c=0)
    """   
    jednadzbe=[]
    for i in range(lijevo.nrows()):
        for j in range(lijevo.ncols()):
            jednadzbe.append(lijevo[i,j]==desno[i,j])
    rj=solve(jednadzbe,*varijable,solution_dict=rjecnik)
    return rj


#### 1. poglavlje ####

def ProjekcijaVektorVektor(vek1,vek2,full=False):
    """Daje ortogonalnu projekciju vektora vek1 na vektor vek2.
       Vektori mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje rezultat pomocu simplify_full ako je to moguce.
       
       -- PRIMJERI --
       var('a b')
       ProjekcijaVektorVektor((1,2,3),vector((0,4,-5)))
       ProjekcijaVektorVektor([1,b,3],(a,-2,1))
       ProjekcijaVektorVektor([1,b,3],(a,-2,1),True)"""
    v1=vector(vek1)
    v2=vector(vek2)
    d2=v2.norm()
    if d2==0: return "Error: projiciranje na nulvektor"
    rez=v1.dot_product(v2)/(d2^2)*v2
    if full:
        try:
            return rez.simplify_full()
        except:
            return rez
    else:
        return rez


def ProjekcijaVektorRavnina(vek,nor,full=False):
    """Daje ortogonalnu projekciju vektora vek na ravninu koja ima za normalu vektor nor.
       Vektori mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje rezultat pomocu simplify_full ako je to moguce.
       
       -- PRIMJERI --
       var('a b')
       ProjekcijaVektorRavnina((1,2,3),vector((0,4,-5)))
       ProjekcijaVektorRavnina([1,b,3],(a,-2,1))
       ProjekcijaVektorRavnina([1,b,3],(a,-2,1),True)"""
    v1=vector(vek)
    v2=vector(nor)
    d2=v2.norm()
    if d2==0: return "Error: normala ravnine je nulvektor"
    rez=v1-v1.dot_product(v2)/(d2^2)*v2
    if full:
        try:
            return rez.simplify_full()
        except:
            return rez
    else:
        return rez


def KutVektora(vek1,vek2,full=False,izlaz='DMS'):
    """Daje kut izmedju vektora vek1 i vek2.
       Vektori mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za kut.
       Opcija full=True pojednostavljuje izraz za kut pomocu simplify_full.
       Opcija izlaz='DMS' (default) daje kut u stupnjevima, minutama, sekundama i stotinkama sekundi.
       Opcija izlaz=None daje egzaktni izraz za kut u radijanima.
       
       -- PRIMJERI --
       var('a')
       KutVektora((4,3,1),(4,1,1))
       KutVektora((4,3,1),(4,1,1),izlaz=None)
       KutVektora((4,3,1),(a,1,1),izlaz=None)"""
    v1=vector(vek1)
    v2=vector(vek2)
    if v1==vector((0,0,0)) or v2==vector((0,0,0)): return "Error: barem jedan od vektora je nulvektor"
    kut=v1.dot_product(v2)/(v1.norm()*v2.norm())
    if full:
        try:
            kut=kut.simplify_full()
        except:
            pass
    if izlaz=='DMS':
        return DMS(acos(kut),'radijan')
    elif izlaz==None:
        return acos(kut)
    else:
        return "Error: izlaz mora biti 'DMS' ili None"


def SP(vek1,vek2,full=False):
    """Daje skalarni produkt vektora vek1 i vek2.
       Vektori mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje po potrebi rezultat pomocu 
       simplify_full ukoliko je to moguce.
       
       -- PRIMJERI --
       SP((1,2,4),(3,-1,0))
       var('a')
       SP((a,2,1),(-1,0,a^2))"""
    v1=vector(vek1)
    v2=vector(vek2)
    rez=v1.dot_product(v2)
    if full:
        try:
            rez=rez.simplify_full()
        except:
            pass
    return rez


def VP(vek1,vek2,full=False):
    """Daje vektorski produkt vektora vek1 i vek2.
       Vektori mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje po potrebi rezultat pomocu 
       simplify_full ukoliko je to moguce.
       
       -- PRIMJERI --
       VP((1,2,4),(3,-1,0))
       var('a')
       VP((a,2,1),(-1,0,a^2))"""
    v1=vector(vek1)
    v2=vector(vek2)
    rez=v1.cross_product(v2)
    if full:
        try:
            rez=rez.simplify_full()
        except:
            pass
    return rez


def MP(vek1,vek2,vek3,full=False):
    """Daje mjesoviti produkt vektora vek1, vek2 i vek3.
       Vektori mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje po potrebi rezultat pomocu 
       simplify_full ukoliko je to moguce.
       
       -- PRIMJERI --
       MP((1,2,4),(3,-1,0),(1,1,-2))
       var('a')
       MP((a,2,1),(-1,0,a^2),(1,1,-2))"""
    v1=vector(vek1)
    v2=vector(vek2)
    v3=vector(vek3)
    rez=v1.cross_product(v2).dot_product(v3)
    if full:
        try:
            rez=rez.simplify_full()
        except:
            pass
    return rez


#### 2. poglavlje ####

def dTocka(T1,T2,full=False):
    """Daje udaljenost tocaka T1 i T2. Tocke mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za udaljenost.
       Opcija full=True pojednostavljuje izraz za udaljenost pomocu simplify_full.
       
       -- PRIMJERI --
       var('a')
       T1=(0,-3,2)
       T2=[4,0,-5]
       T3=vector((2,1,-8))
       dTocka(T1,T2)
       dTocka(T3,(a,7,2))
       dTocka(T3,(a,7,2),True)"""
    v1=vector(T1)
    v2=vector(T2)
    rez=(v1-v2).norm()
    if full:
        try:
            return rez.simplify_full()
        except:
            return rez
    else:
        return rez


def Omjer(T1,T2,lam=-1,full=False):
    """Daje tocku na pravcu T1T2 koja duzinu T1T2 dijeli u omjeru lam.
       Po defaultu je lam=-1, sto daje poloviste duzine T1T2.
       Ako je lam=1, dobivamo beskonacno daleku tocku na pravcu T1T2.
       Tocke mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje rezultat pomocu simplify_full ako je to moguce.
       
       -- PRIMJERI --
       var('a b')
       T1=(0,1,-5)
       T2=vector([4,3,1])
       T3=vector((sqrt(a^2),0,a+b))+vector((1,1,b^2))
       Omjer(T1,T2,lam=2/3)
       Omjer(T3,T1,lam=-4)
       Omjer(T3,T1,lam=-4,full=True)"""
    if lam==1: return Infinity
    v1=vector(T1)
    v2=vector(T2)
    rez=1/(1-lam)*(v1-lam*v2)
    if full:
        try:
            return rez.simplify_full()
        except:
            return rez
    else:
        return rez


def Ravnina3Tocke(T1,T2,T3,skrati=True,ispis=1,full=False):
    """Daje jednadzbu ravnine zadanu s tri nekolinearne tocke T1, T2 i T3.
       Tocke mogu biti tipa list, tuple, vector.
       Opcija ispis=1 (default) daje ravninu u obliku Ax+By+Cz+D=0.
       Opcija ispis=0 daje ravninu u obliku [A,B,C,D].
       Opcija skrati=True (default) u slucaju racionalnih brojeva A,B,C,D vraca oblik
       u kojemu su A,B,C,D do kraja skraceni.
       Opcija skrati=False ne skracuje brojeve A,B,C,D.
       Opcija full=True u slucaju simbolickih izraza pojednostavljuje izraze
       za A,B,C,D pomocu simplify_full ukoliko je to moguce.
       Po defaultu je stavljeno full=False.
       
       -- PRIMJERI --
       T1=(3,2,1)
       T2=[4,0,-1)
       T3=vector((2,2,-5))
       Ravnina3Tocke(T1,T2,T3)"""        
    v1=vector(T1)
    v2=vector(T2)
    v3=vector(T3)
    v12=v2-v1
    v13=v3-v1
    normala=v12.cross_product(v13)
    if normala==vector((0,0,0)): return "Error: tocke su kolinearne"
    D=-normala.dot_product(v1)
    ravnina=normala.list()+[D]
    if full:
        try:
            ravnina=vector(ravnina).simplify_full().list()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(ravnina)
            ravnina=list(map(lambda x: br*x, ravnina))
        except:
            pass
    if ispis==0:
        return ravnina
    if ispis==1:
        var('x y z')
        return vector(ravnina[:-1]).dot_product(vector((x,y,z)))+ravnina[-1]==0


def Ravnina2Vektora(T1,vek1,vek2,skrati=True,ispis=1,full=False):
    """Daje jednadzbu ravnine zadanu s tockom T1 i dva nekolinearna vektora vek1 i vek2.
       Tocka i vektori mogu biti tipa list, tuple, vector.
       Opcija ispis=1 (default) daje ravninu u obliku Ax+By+Cz+D=0.
       Opcija ispis=0 daje ravninu u obliku [A,B,C,D].
       Opcija skrati=True (default) u slucaju racionalnih brojeva A,B,C,D vraca oblik
       u kojemu su A,B,C,D do kraja skraceni.
       Opcija skrati=False ne skracuje brojeve A,B,C,D.
       Opcija full=True u slucaju simbolickih izraza pojednostavljuje izraze
       za A,B,C,D pomocu simplify_full ukoliko je to moguce.
       Po defaultu je stavljeno full=False.
       
       -- PRIMJERI --
       Ravnina2Vektora((1,2,-4),(2,1,2),(3,0,1))
       Ravnina2Vektora((1,2,-4),(2,1,2),[3,0,1],ispis=0)"""  
    v0=vector(T1)
    v1=vector(vek1)
    v2=vector(vek2)
    normala=v1.cross_product(v2)
    if normala==vector((0,0,0)): return "Error: vektori su kolinearni"
    D=-normala.dot_product(v0)
    ravnina=normala.list()+[D]
    if full:
        try:
            ravnina=vector(ravnina).simplify_full().list()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(ravnina)
            ravnina=list(map(lambda x: br*x, ravnina))
        except:
            pass
    if ispis==0:
        return ravnina
    if ispis==1:
        var('x y z')
        return vector(ravnina[:-1]).dot_product(vector((x,y,z)))+ravnina[-1]==0


def RavninaNormala(T1,nor,skrati=True,ispis=1,full=False):
    """Daje jednadzbu ravnine zadanu s tockom T1 i normalom nor.
       Tocka i normala mogu biti tipa list, tuple, vector.
       Opcija ispis=1 (default) daje ravninu u obliku Ax+By+Cz+D=0.
       Opcija ispis=0 daje ravninu u obliku [A,B,C,D].
       Opcija skrati=True (default) u slucaju racionalnih brojeva A,B,C,D vraca oblik
       u kojemu su A,B,C,D do kraja skraceni.
       Opcija skrati=False ne skracuje brojeve A,B,C,D.
       Opcija full=True u slucaju simbolickih izraza pojednostavljuje izraze
       za A,B,C,D pomocu simplify_full ukoliko je to moguce.
       Po defaultu je stavljeno full=False.
       
       -- PRIMJERI --
       RavninaNormala((1,2,-4),(2,1,2))
       RavninaNormala((1,2,-4),(2,1,2),ispis=0)"""
    v1=vector(T1)
    v2=vector(nor)
    if v2==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    D=-v2.dot_product(v1)
    ravnina=v2.list()+[D]
    if full:
        try:
            ravnina=vector(ravnina).simplify_full().list()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(ravnina)
            ravnina=list(map(lambda x: br*x, ravnina))
        except:
            pass
    if ispis==0:
        return ravnina
    if ispis==1:
        var('x y z')
        return vector(ravnina[:-1]).dot_product(vector((x,y,z)))+ravnina[-1]==0


def dTockaRavnina(T1,rav,full=False):
    """Daje udaljenost tocke T1 od ravnine rav. 
       Tocka i ravnina mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za udaljenost.
       Opcija full=True pojednostavljuje izraz za udaljenost pomocu simplify_full.
       
       -- PRIMJERI --
       dTockaRavnina((3,4,1),[2,1,0,-5])"""
    v1=vector(T1)
    normala=vector(rav[:-1])
    if normala==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    rez=abs(normala.dot_product(v1)+rav[-1])/normala.norm()
    if full:
        try:
            return rez.simplify_full()
        except:
            return rez
    else:
        return rez


def KutRavnina(nor1,nor2,full=False,izlaz='DMS'):
    """Daje kut izmedju ravnina s normalama nor1 i nor2.
       Normale mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za kut.
       Opcija full=True pojednostavljuje izraz za kut pomocu simplify_full.
       Opcija izlaz='DMS' (default) daje kut u stupnjevima, minutama, sekundama i stotinkama sekundi.
       Opcija izlaz=None daje egzaktni izraz za kut u radijanima.
       
       -- PRIMJERI --
       var('a')
       KutRavnina((4,3,1),(4,1,1))
       KutRavnina((4,3,1),(4,1,1),izlaz=None)
       KutRavnina((4,3,1),(a,1,1),izlaz=None)"""
    v1=vector(nor1)
    v2=vector(nor2)
    if v1==vector((0,0,0)) or v2==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    kut=abs(v1.dot_product(v2))/(v1.norm()*v2.norm())
    if full:
        try:
            kut=kut.simplify_full()
        except:
            pass
    if izlaz=='DMS':
        return DMS(acos(kut),'radijan')
    elif izlaz==None:
        return acos(kut)
    else:
        return "Error: izlaz mora biti 'DMS' ili None"


def RavninaNormalniOblik(rav,ispis=1,full=False):
    """Daje normalni oblik jednadzbe ravnine.
       Ravnina moze biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje rezultat pomocu simplify_full.
       Opcija ispis=1 (default) daje normalni oblik u obliku Ax+By+Cz+D=0.
       Opcija ispis=0 daje normalni oblik u obliku [A,B,C,D].
       
       -- PRIMJERI --
       RavninaNormalniOblik([1,2,3,-2])
       RavninaNormalniOblik([1,2,3,-2],ispis=0)""" 
    vek=vector(rav)
    duljina=vek[:-1].norm()
    if duljina==0: return "Error: normala ravnine je nulvektor"
    if vek[-1]==0: 
        ravnina=(1/duljina)*vek
    else:
        ravnina=-1/(sign(vek[-1])*duljina)*vek
    if full:
        try:
            ravnina=ravnina.simplify_full()
        except:
            pass
    if ispis==0:
        return ravnina.list()
    if ispis==1:
        var('x y z')
        return ravnina[:-1].dot_product(vector((x,y,z)))+ravnina[-1]==0


def Pravac2Tocke(T1,T2,skrati=True,full=False):
    """Daje vektor smjera pravca odredjenog tockama T1 i T2.
       Tocke mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno vektor smjera.
       Opcija full=True pojednostavljuje vektor smjera pomocu simplify_full.
       Opcija skrati=True (default) u slucaju racionalnih brojeva vraca
       vektor smjera (a,b,c) tako da je GCD(a,b,c)=1.
       Opcija skrati=False vraca samo vektor smjera bez dodatnog kracenja komponenti.
       
       -- PRIMJERI --
       var('a')
       Pravac2Tocke((1,2,-3),(3,-2,5))
       Pravac2Tocke((sqrt(a^2),4,-2),(5,8,8),full=True)"""       
    v1=vector(T1)
    v2=vector(T2)
    s=v1-v2
    if s==vector((0,0,0)): return "Error: pravac nije odredjen"
    if full:
        try:
            s=s.simplify_full()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(s.list())
            s=s.apply_map(lambda x: br*x)
        except:
            pass
    return s


def PresjekRavnina(rav1,rav2,tocka=None,skrati=True,full=False):
    """Daje presjek ravnina rav1 i rav2. 
       Ravnine mogu biti tipa list, tuple, vector. Ravnine se unose tako da
       se unose redom njihovi koeficijenti iz opceg oblika.
       Ukoliko se ravnine podudaraju ili su paralelne, vraca se odgovarajuca poruka.
       Ukoliko se ravnine sijeku po pravcu, vraca se rjecnik s vektorom smjera i jednom
       tockom na tom pravcu.
       
       Opcija full=False (default) ne pojednostavljuje dodatno vektor smjera.
       Opcija full=True pojednostavljuje vektor smjera pomocu simplify_full.
       Opcija skrati=True (default) u slucaju racionalnih brojeva vraca
       vektor smjera (a,b,c) tako da je GCD(a,b,c)=1.
       Opcija skrati=False vraca samo vektor smjera bez dodatnog kracenja komponenti.
       
       Ukoliko je tocka=None, tada se pokusava pronaci neka tocka na pravcu tako da
       se prvo proba s x=0, pa ako to ne uspije proba se s y=0, pa ako to ne uspije
       na kraju se stavi z=0.
       
       Ako stavimo tocka={1:5}, tada to znaci da zelimo da je x=5.
       Ako stavimo tocka={2:-1}, tada to znaci da zelimo da je y=-1.
       Ako stavimo tocka={3:0}, tada to znaci da zelimo da je z=0.
       Dakle, opcija tocka je rjecnik sa samo jednim kljucem pri cemu kljuc moze
       biti 1,2 ili 3. Ovi brojevi redom oznacavaju koordinate x,y i z.
       Vrijednost kljuca je zapravo vrijednost koju zelimo da ima odgovarajuca koordinata,
       kako je gore pokazano na konkretnim primjerima.
       Ukoliko se ne uspije pronaci tocka na presjecnom pravcu koja bi zadovoljavala nas
       postavljeni zahtjev, tada algoritam automatski prelazi na opciju tocka=None.
       
       -- PRIMJERI --
       PresjekRavnina([3,1,-17,0],[2,3,-16,-7])
       PresjekRavnina([3,1,-17,0],[2,3,-16,-7],tocka={3:1})
       PresjekRavnina([1,0,0,-3],[0,1,0,-2],tocka={1:0})"""       
    n1=vector(rav1[:-1])
    n2=vector(rav2[:-1])
    if n1==vector((0,0,0)): return "Error: normala prve ravnine je nulvektor"
    if n2==vector((0,0,0)): return "Error: normala druge ravnine je nulvektor"
    M1=matrix([rav1[:-1],rav2[:-1]])
    M2=matrix([rav1,rav2])
    if M2.rank()==1: return "Ravnine su jednake."
    if M1.rank()==1: return "Ravnine su paralelne."
    s=n1.cross_product(n2)
    if full:
        try:
            s=s.simplify_full()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(s.list())
            s=s.apply_map(lambda x: br*x)
        except:
            pass
    var('x y z')
    nepoznanice=[x,y,z]
    sustav=[n1.dot_product(vector((x,y,z)))+rav1[-1]==0, n2.dot_product(vector((x,y,z)))+rav2[-1]==0]
    if tocka != None:
        xyz=nepoznanice[list(tocka.keys())[0]-1]
        nepoznanice2=list(filter(lambda q: q!=xyz, nepoznanice))
        vrijednost=list(tocka.values())[0]
        try:
           rjesenje=solve(list(map(lambda q: q.substitute({xyz:vrijednost}), sustav)),*nepoznanice2,solution_dict=True)[0]
           rjesenje[xyz]=vrijednost
           tocka_pravac=vector([rjesenje[q] for q in nepoznanice])
           pravac={'tocka':tocka_pravac, 'vektor_smjera':s}
           return pravac
        except:
            print("Nije uspjelo sa zadanom pocetnom vrijednosti varijable.")
            for nep in nepoznanice:
                try:
                    nepoznanice2=list(filter(lambda q: q!=nep, nepoznanice))
                    rjesenje=solve(list(map(lambda q: q.substitute({nep:0}), sustav)),*nepoznanice2,solution_dict=True)[0]
                    rjesenje[nep]=0
                    tocka_pravac=vector([rjesenje[q] for q in nepoznanice])
                    pravac={'tocka':tocka_pravac, 'vektor_smjera':s}
                    return pravac
                except:
                    pass        
    else:
        for nep in nepoznanice:
            try:
                nepoznanice2=list(filter(lambda q: q!=nep, nepoznanice))
                rjesenje=solve(list(map(lambda q: q.substitute({nep:0}), sustav)),*nepoznanice2,solution_dict=True)[0]
                rjesenje[nep]=0
                tocka_pravac=vector([rjesenje[q] for q in nepoznanice])
                pravac={'tocka':tocka_pravac, 'vektor_smjera':s}
                return pravac
            except:
                pass


def KutPravaca(vek1,vek2,full=False,izlaz='DMS'):
    """Daje kut izmedju pravaca s vektorima smjerova vek1 i vek2.
       Vektori smjerova mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za kut.
       Opcija full=True pojednostavljuje izraz za kut pomocu simplify_full.
       Opcija izlaz='DMS' (default) daje kut u stupnjevima, minutama, sekundama i stotinkama sekundi.
       Opcija izlaz=None daje egzaktni izraz za kut u radijanima.
       
       -- PRIMJERI --
       var('a')
       KutPravaca((4,3,1),(4,1,1))
       KutPravaca((4,3,1),(4,1,1),izlaz=None)
       KutPravaca((4,3,1),(a,1,1),izlaz=None)"""
    v1=vector(vek1)
    v2=vector(vek2)
    if v1==vector((0,0,0)) or v2==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    kut=abs(v1.dot_product(v2))/(v1.norm()*v2.norm())
    if full:
        try:
            kut=kut.simplify_full()
        except:
            pass
    if izlaz=='DMS':
        return DMS(acos(kut),'radijan')
    elif izlaz==None:
        return acos(kut)
    else:
        return "Error: izlaz mora biti 'DMS' ili None"


def KutRavninaPravac(nor,vek,full=False,izlaz='DMS'):
    """Daje kut izmedju ravnine s normalom nor i pravca s vektorom smjera vek.
       Normala ravnine i vektor smjera pravca mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za kut.
       Opcija full=True pojednostavljuje izraz za kut pomocu simplify_full.
       Opcija izlaz='DMS' (default) daje kut u stupnjevima, minutama, sekundama i stotinkama sekundi.
       Opcija izlaz=None daje egzaktni izraz za kut u radijanima.
       
       -- PRIMJERI --
       var('a')
       KutRavninaPravac((4,3,1),(4,1,1))
       KutRavninaPravac((4,3,1),(4,1,1),izlaz=None)
       KutRavninaPravac((4,3,1),(a,1,1),izlaz=None)"""
    v1=vector(nor)
    v2=vector(vek)
    if v1==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    if v2==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    kut=abs(v1.dot_product(v2))/(v1.norm()*v2.norm())
    if full:
        try:
            kut=kut.simplify_full()
        except:
            pass
    if izlaz=='DMS':
        return DMS(asin(kut),'radijan')
    elif izlaz==None:
        return asin(kut)
    else:
        return "Error: izlaz mora biti 'DMS' ili None"


def dTockaPravac(T1,vek,T0,full=False):
    """Daje udaljenost tocke T1 od pravca koji prolazi tockom T0 i ima vektor smjera vek. 
       Tocke i vektor smjera mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno izraz za udaljenost.
       Opcija full=True pojednostavljuje izraz za udaljenost pomocu simplify_full.
       
       -- PRIMJERI --
       var('a')
       dTockaPravac((6,1,-5),(3,1,-1),(2,-4,2))
       dTockaPravac((a,1,-5),(3,1,-1),(2,-4,2),full=True)"""
    r1=vector(T1)
    s=vector(vek)
    r0=vector(T0)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    brojnik=(r1-r0).cross_product(s).norm()
    nazivnik=s.norm()
    rez=brojnik/nazivnik
    if full:
        try:
            rez=rez.simplify_full()
        except:
            pass
    return rez


def PresjekPravaca(vek1,T1,vek2,T2,full=False):
    """Daje presjek pravca koji prolazi tockom T1 i ima vektor smjera vek1
       s pravcem koji prolazi tockom T2 i ima vektor smjera vek2.
       Tocke i vektori smjerova mogu biti tipa list, tuple, vector.
       Ako su pravci paralelni ili mimosmjerni, vraca se odgovarajuca poruka.
       Opcija full=False (default) ne pojednostavljuje dodatno konacni rezultat.
       Opcija full=True pojednostavljuje konacni rezultat pomocu simplify_full.
       
       Ukoliko se pravci sijeku, na izlazu se vraca lista od tri elementa.
       Prvi element liste je tocka presjeka, drugi element je vrijednost parametra
       u za koji se dobiva tocka presjeka na prvom pravcu, a treci element je vrijednost
       parametra v za koji se dobiva tocka presjeka na drugom pravcu.
       
       -- PRIMJERI --
       PresjekPravaca((8,-3,1),(1,-2,-1),(2,1,3),(6,-3,1))"""
    v1=vector(vek1)
    v2=vector(vek2)
    if v1==vector((0,0,0)) or v2==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    r1=vector(T1)
    r2=vector(T2)
    s=v1.cross_product(v2)
    r12=v1.cross_product(r1-r2)
    if s==vector((0,0,0)) and r12==vector((0,0,0)): return "Pravci su jednaki"
    if s==vector((0,0,0)): return "Pravci su paralelni"
    br=matrix([r2-r1,v1,v2]).determinant()
    if br!=0: return "Pravci su mimosmjerni"
    var('x y z u v')
    sustav=[x==r1[0]+v1[0]*u,y==r1[1]+v1[1]*u,z==r1[2]+v1[2]*u,x==r2[0]+v2[0]*v,y==r2[1]+v2[1]*v,z==r2[2]+v2[2]*v]
    rez=solve(sustav,x,y,z,u,v,solution_dict=True)[0]
    tocka=vector((rez[x],rez[y],rez[z]))
    rez_u=rez[u]
    rez_v=rez[v]
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
        try:
            rez_u=rez_u.simplify_full()
        except:
            pass
        try:
            rez_v=rez_v.simplify_full()
        except:
            pass
    rez=[tocka,u==rez_u,v==rez_v]
    return rez


def PresjekPravacRavnina(vek,T0,rav,full=False):
    """Daje presjek pravca koji prolazi tockom T0 i ima vektor smjera vek s ravninom 
       rav koja se zadaje preko svojih koeficijenata iz opceg oblika jednadzbe.
       Ukoliko su pravac i ravnina paralelni, vraca se odgovarajuca poruka.
       Ulazni podaci za pravac i ravninu mogu biti tipa list, tuple, vector.
       Opcija full=False (default) ne pojednostavljuje dodatno konacni rezultat.
       Opcija full=True pojednostavljuje konacni rezultat pomocu simplify_full.
       
       Ukoliko se pravac i ravnina sijeku, na izlazu se vraca lista od dva elementa.
       Prvi element je tocka presjeka, a drugi element je vrijednost parametra t za
       koji se tocka presjeka dobiva na zadanom pravcu.
       
       -- PRIMJERI --
       PresjekPravacRavnina((3,1,-1),(2,-4,2),[3,1,-1,-24])"""
    s=vector(vek)
    nor=vector(rav[:-1])
    r0=vector(T0)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    if nor==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    sp=s.dot_product(nor)
    br=nor.dot_product(r0)+rav[-1]
    if sp==0 and br==0: return "Pravac lezi u ravnini"
    if sp==0: return "Pravac i ravnina su paralelni"
    var('x y z t')
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r0[0]+s[0]*t,y==r0[1]+s[1]*t,z==r0[2]+s[2]*t]
    rez=solve(sustav,x,y,z,t,solution_dict=True)[0]
    tocka=vector((rez[x],rez[y],rez[z]))
    rez_t=rez[t]
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
        try:
            rez_t=rez_t.simplify_full()
        except:
            pass
    rez=[tocka,t==rez_t]
    return rez


def ZajednickaNormala(vek1,T1,vek2,T2,izlaz=1,skrati=True,full=False):
    """Daje zajednicku normalu mimosmjernih pravaca. 
       Prvi pravac ima vektor smjera vek1 i prolazi tockom T1, a drugi
       pravac ima vektor smjera vek2 i prolazi tockom T2.
       Tocke i vektori mogu biti tipa list, tuple, vector.
       
       Ako su pravci paralelni, vraca se error.
       Ako se pravci sijeku, na izlazu se vraca zajednicka normala kao
       rjecnik s vektorom smjera i presjecnom tockom zadana dva pravca.
       Ispisuje se takodjer poruka da se pravci sijeku.
       
       Ako su pravci mimosmjerni, tada se na izlazu zajednicka normala
       daje u jednom od dva oblika.
       Ako je izlaz=1, tada se vraca rjecnik sa vektorom smjera normale, 
       presjecnim tockama tocka1 i tocka2 normale redom sa zadanim pravcima i 
       sa vrijednostima parametara u i v za koje se te tocke postizu redom 
       na pocetnim pravcima.
       Ako je izlaz=2, tada se zajednicka normala daje kao presjek dvije ravnine,
       a ravnine su dane u listi preko svojih opcih jednadzbi.
       Ako je izlaz=3, tada se zajednicka normala daje kao presjek dvije ravnine,
       a ravnine su dane u listi kao liste koeficijenata iz svojih opcih oblika 
       jednadzbi.
              
       Opcija full=False (default) ne pojednostavljuje dodatno rezultate.
       Opcija full=True pojednostavljuje dodatno rezultate pomocu full_simplify
       ako je to moguce. To moze biti od koristi ako na ulazu ima simbolickih izraza.
       
       Opcija skrati=True (default) skracuje komponente vektora smjera zajednicke normale
       i koeficijente ravnina u slucaju da je izlaz=2, ukoliko je to moguce.
       Opcija skrati=False nece nista skracivati.
       
       -- PRIMJERI --
       var('a')
       ZajednickaNormala((1,1,2),(0,0,1),(1,2,1),(1,0,0))
       ZajednickaNormala((1,1,2),(0,0,1),(1,2,1),(1,0,0),izlaz=2)
       ZajednickaNormala((8,-3,1),(1,-2,-1),(2,1,3),(6,-3,1))
       ZajednickaNormala((a,-3,1),(1,-2,-1),(2,1,3),(6,-3,1),izlaz=2,full=True)
       ZajednickaNormala((a,-3,1),(1,-2,-1),(2,1,3),(6,-3,1),izlaz=3,full=True)"""
    v1=vector(vek1)
    v2=vector(vek2)
    r1=vector(T1)
    r2=vector(T2)
    if v1==vector((0,0,0)) or v2==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    s=v1.cross_product(v2)
    if s==vector((0,0,0)): return "Error: pravci su paralelni"
    if full:
        try:
            s=s.simplify_full()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(s.list())
            s=s.apply_map(lambda x: br*x)
        except:
            pass
    br=matrix([r2-r1,v1,v2]).determinant()
    if br==0:
        var('x y z u v')
        sustav=[x==r1[0]+v1[0]*u,y==r1[1]+v1[1]*u,z==r1[2]+v1[2]*u,x==r2[0]+v2[0]*v,y==r2[1]+v2[1]*v,z==r2[2]+v2[2]*v]
        rez=solve(sustav,x,y,z,u,v,solution_dict=True)[0]
        tocka=vector((rez[x],rez[y],rez[z])) 
        if full:
            try:
                tocka=tocka.simplify_full()
            except:
                pass
        pravac={'tocka':tocka,'vektor_smjera':s}
        print("Pravci se sijeku")
        return pravac
    if izlaz==1:
        var('u v k')
        vektor1=(r2+v*v2)-(r1+u*v1)
        vektor2=k*s
        rez=solve([vektor1[0]==vektor2[0],vektor1[1]==vektor2[1],vektor1[2]==vektor2[2]],u,v,k,solution_dict=True)[0]
        tocka1=(r1+u*v1).subs(u==rez[u])
        tocka2=(r2+v*v2).subs(v==rez[v])
        if full:
            try:
                tocka1=tocka1.simplify_full()
            except:
                pass
            try:
                tocka2=tocka2.simplify_full()
            except:
                pass
        return {'tocka1':tocka1,'tocka2':tocka2,'vektor_smjera':s,u:rez[u],v:rez[v]}
    elif izlaz==2 or izlaz==3:
        var('x y z')
        r=vector((x,y,z))
        nor1=v1.cross_product(s)
        D1=-nor1.dot_product(r1)
        ravnina1=nor1.list()+[D1]
        nor2=v2.cross_product(s)
        D2=-nor2.dot_product(r2)
        ravnina2=nor2.list()+[D2]
        if full:
            try:
                ravnina1=vector(ravnina1).simplify_full().list()
            except:
                pass
            try:
                ravnina2=vector(ravnina2).simplify_full().list()
            except:
                pass
        if skrati:
            try:
                br=1/gcd(ravnina1)
                ravnina1=list(map(lambda x: br*x, ravnina1))
            except:
                pass
            try:
                br=1/gcd(ravnina2)
                ravnina2=list(map(lambda x: br*x, ravnina2))
            except:
                pass
        rez1=vector(ravnina1[:-1]).dot_product(r)+ravnina1[-1]
        rez2=vector(ravnina2[:-1]).dot_product(r)+ravnina2[-1]
        if izlaz==2:
            return [rez1==0,rez2==0]
        else:
            return [ravnina1,ravnina2]
    else:
        return "Error: 'izlaz' mora biti jednak 1 ili 2 ili 3"


def UdaljenostPravaca(vek1,T1,vek2,T2,full=False):
    """Daje udaljenost pravaca od kojih prvi prolazi tockom T1 i ima vektor smjera vek1,
       a drugi prolazi tockom T2 i ima vektor smjera vek2.
       Tocke i vektori smjerova mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje rezultat pomocu simplify_full ukoliko je to moguce.
       
       -- PRIMJERI --
       UdaljenostPravaca((1,1,2),(0,0,1),(1,2,1),(1,0,0))
       var('a')
       UdaljenostPravaca((2,2,-1),(a,-1,3),(2,2,-1),(0,1,-4),True)"""
    v1=vector(vek1)
    v2=vector(vek2)
    r1=vector(T1)
    r2=vector(T2)
    if v1==vector((0,0,0)) or v2==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    s=v1.cross_product(v2)
    if s==vector((0,0,0)):
        return dTockaPravac(T1,vek2,T2,full)
    else:
        br=abs(matrix([r2-r1,v1,v2]).determinant())
        br2=br/s.norm()
        if full:
            try:
                br2=br2.simplify_full()
            except:
                pass
        return br2


def ProjekcijaTockaRavnina(T0,rav,full=False):
    """Daje ortogonalnu projekciju tocke T0 na ravninu rav.
       Ravnina se zadaje preko svojih koeficijenata iz opceg 
       oblika jednadzbe.
       Tocka i ravnina mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce.
       
       Na izlazu se zapravo daje lista od dva elementa. Prvi element je ortogonalna
       projekcija tocke na ravninu, a drugi element je vrijednost parametra za koji
       se ta tocka dobiva kao presjek zadane ravnine i pravca koji prolazi tockom T0
       i ima za vektor smjera normalu zadane ravnine.
       
       -- PRIMJERI --
       ProjekcijaTockaRavnina((1,2,3),[2,1,-1,-13])
       var('a')
       ProjekcijaTockaRavnina((1,a,3),[2,1,-1,-13])"""   
    r0=vector(T0)
    nor=vector(rav[:-1])
    if nor==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    var('x y z t')
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r0[0]+nor[0]*t,y==r0[1]+nor[1]*t,z==r0[2]+nor[2]*t]
    rez=solve(sustav,x,y,z,t,solution_dict=True)[0]
    tocka=vector((rez[x],rez[y],rez[z]))
    rez_t=rez[t]
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
        try:
            rez_t=rez_t.simplify_full()
        except:
            pass
    rez=[tocka,t==rez_t]
    return rez


def SimetricnaTockaRavnina(T0,rav,full=False):
    """Daje simetricnu tocku tocke T0 s obzirom na ravninu rav.
       Ravnina se zadaje preko svojih koeficijenata iz opceg 
       oblika jednadzbe.
       Tocka i ravnina mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce.
       
       Na izlazu se zapravo daje lista od dva elementa. Prvi element je simetricna
       tocka, a drugi element je vrijednost parametra za koji se ta tocka dobiva na 
       pravcu koji prolazi tockom T0 i ima za vektor smjera normalu zadane ravnine.
       
       -- PRIMJERI --
       SimetricnaTockaRavnina((1,2,3),[2,1,-1,-13])
       var('a')
       SimetricnaTockaRavnina((1,a,3),[2,1,-1,-13])"""   
    r0=vector(T0)
    nor=vector(rav[:-1])
    if nor==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    var('x y z t')
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r0[0]+nor[0]*t,y==r0[1]+nor[1]*t,z==r0[2]+nor[2]*t]
    rez_t=2*solve(sustav,x,y,z,t,solution_dict=True)[0][t]
    tocka=vector((r0[0]+nor[0]*rez_t,r0[1]+nor[1]*rez_t,r0[2]+nor[2]*rez_t))
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
        try:
            rez_t=rez_t.simplify_full()
        except:
            pass
    rez=[tocka,t==rez_t]
    return rez


def ProjekcijaPravacRavnina(vek,T1,rav,skrati=True,full=False):
    """Daje ortogonalnu projekciju pravca koji prolazi tockom T1 i 
       ima vektor smjera vek s obzirom na ravninu rav.
       Ravnina se zadaje preko svojih koeficijenata iz opceg 
       oblika jednadzbe.
       Tocka, vektor smjera i ravnina mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce.
       
       Opcija skrati=True (default) u slucaju racionalnih brojeva vraca
       vektor smjera (a,b,c) ortogonalne projekcije tako da je GCD(a,b,c)=1.
       Opcija skrati=False vraca samo vektor smjera bez dodatnog kracenja komponenti.
       
       Na izlazu se vraca pravac kao rjecnik s tockom i vektorom smjera ortogonalne
       projekcije.
       
       -- PRIMJERI --
       ProjekcijaPravacRavnina((2,1,-1),(-1,0,0),[1,2,-2,-5])
       var('a')
       ProjekcijaPravacRavnina((a,1,-1),(-1,0,0),[1,2,-2,-5],full=True)"""   
    s=vector(vek)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    nor=vector(rav[:-1])
    if nor==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    r1=vector(T1)
    var('x y z t')
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r1[0]+nor[0]*t,y==r1[1]+nor[1]*t,z==r1[2]+nor[2]*t]
    rez_t=solve(sustav,x,y,z,t,solution_dict=True)[0][t]
    tocka=vector((r1[0]+nor[0]*rez_t,r1[1]+nor[1]*rez_t,r1[2]+nor[2]*rez_t))
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
    if nor.dot_product(s)==0: return {'tocka':tocka,'vektor_smjera':s}
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r1[0]+s[0]*t,y==r1[1]+s[1]*t,z==r1[2]+s[2]*t]
    rez_t=solve(sustav,x,y,z,t,solution_dict=True)[0][t]
    presjek=vector((r1[0]+s[0]*rez_t,r1[1]+s[1]*rez_t,r1[2]+s[2]*rez_t))
    smjer=tocka-presjek
    if full:
        try:
            presjek=presjek.simplify_full()
        except:
            pass
        try:
            smjer=smjer.simplify_full()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(smjer.list())
            smjer=smjer.apply_map(lambda x: br*x)
        except:
            pass
    return {'tocka':presjek,'vektor_smjera':smjer}


def SimetricniPravacRavnina(vek,T1,rav,skrati=True,full=False):
    """Daje simetricni pravac pravca koji prolazi tockom T1 i 
       ima vektor smjera vek s obzirom na ravninu rav.
       Ravnina se zadaje preko svojih koeficijenata iz opceg 
       oblika jednadzbe.
       Tocka, vektor smjera i ravnina mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce.
       
       Opcija skrati=True (default) u slucaju racionalnih brojeva vraca
       vektor smjera (a,b,c) simetricnog pravca tako da je GCD(a,b,c)=1.
       Opcija skrati=False vraca samo vektor smjera bez dodatnog kracenja komponenti.
       
       Na izlazu se vraca pravac kao rjecnik s tockom i vektorom smjera simetricnog
       pravca.
       
       -- PRIMJERI --
       SimetricniPravacRavnina((2,1,-1),(-1,0,0),[1,2,-2,-5])
       var('a')
       SimetricniPravacRavnina((a,1,-1),(-1,0,0),[1,2,-2,-5],full=True)"""   
    s=vector(vek)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    nor=vector(rav[:-1])
    if nor==vector((0,0,0)): return "Error: normala ravnine je nulvektor"
    r1=vector(T1)
    var('x y z t')
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r1[0]+nor[0]*t,y==r1[1]+nor[1]*t,z==r1[2]+nor[2]*t]
    rez_t=2*solve(sustav,x,y,z,t,solution_dict=True)[0][t]
    tocka=vector((r1[0]+nor[0]*rez_t,r1[1]+nor[1]*rez_t,r1[2]+nor[2]*rez_t))
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
    if nor.dot_product(s)==0: return {'tocka':tocka,'vektor_smjera':s}
    sustav=[rav[0]*x+rav[1]*y+rav[2]*z+rav[3]==0,x==r1[0]+s[0]*t,y==r1[1]+s[1]*t,z==r1[2]+s[2]*t]
    rez_t=solve(sustav,x,y,z,t,solution_dict=True)[0][t]
    presjek=vector((r1[0]+s[0]*rez_t,r1[1]+s[1]*rez_t,r1[2]+s[2]*rez_t))
    smjer=tocka-presjek
    if full:
        try:
            presjek=presjek.simplify_full()
        except:
            pass
        try:
            smjer=smjer.simplify_full()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(smjer.list())
            smjer=smjer.apply_map(lambda x: br*x)
        except:
            pass
    return {'tocka':presjek,'vektor_smjera':smjer}


def ProjekcijaTockaPravac(T0,vek,T1,full=False):
    """Daje ortogonalnu projekciju tocke T0 na pravac koji
       prolazi tockom T1 i ima vektor smjera vek.
       Tocke i vektor smjera mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce.
       
       Na izlazu se zapravo daje lista od dva elementa. Prvi element je ortogonalna
       projekcija tocke na zadani pravac, a drugi element je vrijednost parametra za koji
       se ta tocka dobiva na zadanom pravcu.   
       
       -- PRIMJERI --
       var('a')
       ProjekcijaTockaPravac((4,4,3),(4,-2,1),(5,-2,1))
       ProjekcijaTockaPravac((4,4,a),(4,-2,1),(5,-2,1))"""
    r0=vector(T0)
    s=vector(vek)
    r1=vector(T1)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    var('x y z t')
    sustav=[s[0]*(x-r0[0])+s[1]*(y-r0[1])+s[2]*(z-r0[2])==0,x==r1[0]+s[0]*t,y==r1[1]+s[1]*t,z==r1[2]+s[2]*t]
    rez=solve(sustav,x,y,z,t,solution_dict=True)[0]
    tocka=vector((rez[x],rez[y],rez[z]))
    rez_t=rez[t]
    if full:
        try:
            tocka=tocka.simplify_full()
        except:
            pass
        try:
            rez_t=rez_t.simplify_full()
        except:
            pass
    rez=[tocka,t==rez_t]
    return rez


def NormalaTockaPravac(T0,vek,T1,skrati=True,full=False):
    """Daje vektor smjera normale iz tocke T0 na pravac koji
       prolazi tockom T1 i ima vektor smjera vek.
       Tocke i vektor smjera mogu biti tipa list, tuple, vector.
       
       Opcija skrati=True (default) u slucaju racionalnih brojeva vraca
       vektor smjera (a,b,c) tako da je GCD(a,b,c)=1.
       Opcija skrati=False vraca samo vektor smjera bez dodatnog kracenja komponenti.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce. 
       
       -- PRIMJERI --
       var('a')
       NormalaTockaPravac((6,1,-5),(3,1,-1),(2,-4,2))
       NormalaTockaPravac((6,1,-5),(3,1,-1),(2,-4,2),skrati=False)
       NormalaTockaPravac((4,4,a),(4,-2,1),(5,-2,1))"""
    r0=vector(T0)
    s=vector(vek)
    r1=vector(T1)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    var('x y z t')
    sustav=[s[0]*(x-r0[0])+s[1]*(y-r0[1])+s[2]*(z-r0[2])==0,x==r1[0]+s[0]*t,y==r1[1]+s[1]*t,z==r1[2]+s[2]*t]
    rez=solve(sustav,x,y,z,t,solution_dict=True)[0]
    tocka=vector((rez[x],rez[y],rez[z]))
    smjer=tocka-r0
    if full:
        try:
            smjer=smjer.simplify_full()
        except:
            pass
    if skrati:
        try:
            br=1/gcd(smjer.list())
            smjer=smjer.apply_map(lambda x: br*x)
        except:
            pass
    return smjer


def SimetricnaTockaPravac(T0,vek,T1,full=False):
    """Daje simetricnu tocku tocke T0 s obzirom na pravac koji
       prolazi tockom T1 i ima vektor smjera vek.
       Tocke i vektor smjera mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       u slucaju opcih brojeva, ukoliko je to moguce.
       
       -- PRIMJERI --
       var('a')
       SimetricnaTockaPravac((4,4,3),(4,-2,1),(5,-2,1))
       SimetricnaTockaPravac((4,4,a),(4,-2,1),(5,-2,1))"""
    r0=vector(T0)
    s=vector(vek)
    r1=vector(T1)
    if s==vector((0,0,0)): return "Error: vektor smjera pravca je nulvektor"
    var('x y z t')
    sustav=[s[0]*(x-r0[0])+s[1]*(y-r0[1])+s[2]*(z-r0[2])==0,x==r1[0]+s[0]*t,y==r1[1]+s[1]*t,z==r1[2]+s[2]*t]
    rez=solve(sustav,x,y,z,t,solution_dict=True)[0]
    tocka=vector((rez[x],rez[y],rez[z]))
    rez=2*tocka-r0
    if full:
        try:
            rez=rez.simplify_full()
        except:
            pass
    return rez


def PovrsinaTrokuta(T1,T2,T3,full=False):
    """Daje povrsinu trokuta kojemu su vrhovi tocke T1,T2 i T3.
       Tocke mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat u slucaju opcih
       parametara pomocu simplify_full naredbe ukoliko je to moguce.
       
       -- PRIMJERI --
       var('a')
       PovrsinaTrokuta((1,0,0),(0,1,2),(-3,0,1))
       PovrsinaTrokuta((a,0,0),(0,1,2),(-3,0,1),True)"""
    r1=vector(T1)
    r2=vector(T2)
    r3=vector(T3)
    v1=r2-r1
    v2=r3-r1
    povrsina=1/2*v1.cross_product(v2).norm()
    if full:
        try:
            povrsina=povrsina.simplify_full()
        except:
            pass
    if povrsina==0: return "Tocke su kolinearne" 
    return povrsina


def OpsegTrokuta(T1,T2,T3,full=False):
    """Daje opseg trokuta s vrhovima u tockama T1,T2 i T3.
       Tocke mogu biti tipa list, tuple, vector.
       
       Opcija full=False (default) ne pojednostavljuje dodatno rezultat.
       Opcija full=True pojednostavljuje dodatno rezultat pomocu simplify_full
       ukoliko je to moguce.
       
       Na izlazu se zapravo daje lista od dva elementa. Prvi element liste je lista 
       redom duljina stranica T2T3, T1T3 i T1T2. Drugi element liste je opseg trokuta.
       
       -- PRIMJERI --
       var('a')
       OpsegTrokuta((1,0,0),(0,1,2),(-3,0,1))
       OpsegTrokuta((a,0,0),(0,1,2),(-3,0,1),True)""" 
    r1=vector(T1)
    r2=vector(T2)
    r3=vector(T3)
    v1=(r2-r3).norm()
    v2=(r1-r3).norm()
    v3=(r1-r2).norm()
    opseg=v1+v2+v3
    w=(r3-r1).cross_product(r2-r1).norm()
    try:
        w=w.simplify_full()
    except:
        pass
    if w==0: return "Tocke su kolinearne" 
    if full:
        try:
            v1=v1.simplify_full()
        except:
            pass
        try:
            v2=v2.simplify_full()
        except:
            pass
        try:
            v3=v3.simplify_full()
        except:
            pass  
        try:
            opseg=opseg.simplify_full()
        except:
            pass
    return [[v1,v2,v3],opseg]


#### 3. poglavlje ####

def matrica_prijelaza(baza1,baza2):
    """Daje matricu prijelaza iz baze baza1 u bazu baza2.
    
    --- PRIMJER ---
     
     B1=[(1,0,0),(1,1,0),(3,4,5)]
     B2=[(3,1,-1),(1,-1,-1),(0,0,2)]
     matrica_prijelaza(B1,B2)
    """
    m1=matrix(baza1).transpose()
    if m1.rank()!=m1.ncols() or m1.rank()!=m1.nrows(): return "Error: prvi skup nije baza"
    m2=matrix(baza2).transpose()
    if m2.rank()!=m2.ncols() or m2.rank()!=m2.nrows(): return "Error: drugi skup nije baza"
    if m1.nrows()!=m2.nrows(): return "Error: baze nisu iz istog vektorskog prostora"
    return m1^-1*m2


#### 4. poglavlje ####

def operator_kan(f):
    """Daje matricni zapis linearnog operatora f u paru kanonskih baza.
    Linearni operator mora biti prethodno definiran kao funkcija ili moze biti
    definiran na ulazu kao lambda funkcija.
    
    --- PRIMJER ---
     
     def f(x,y):
         return (x+y,x-2*y,3*y)
     
     operator_kan(f)
     
     operator_kan(lambda x,y:(x+y,x-2*y,3*y))
    """     
    n=f.__code__.co_argcount
    kan=identity_matrix(n).columns()
    slika=map(lambda t: f(*t), kan)
    return matrix(slika).transpose()


def operator_kanbaza(f,baza1,baza2=None):
    """Daje matricni zapis linearnog operatora f u paru baza (baza1,baza2).
    Ako se baza2 eksplicitno ne navede, podrazumijeva se da je baza2=baza1.
    
    Linearni operator f moze biti:
      * prethodno definiran kao funkcija 
      * definiran na ulazu kao lambda funkcija
      * definiran svojim matricnim zapisom u paru kanonskih baza
    
    --- PRIMJER ---
     
     def f(x,y):
         return (x+y,x-2*y,3*y)
     
     operator_kanbaza(f,[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)])
     
     operator_kanbaza(lambda x,y:(x+y,x-2*y,3*y),[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)])
     
     operator_kanbaza(matrix([[1,1],[1,-2],[0,3]]),[(1,-1),(1,0)],[(0,1,0),(1,2,0),(1,1,1)])
    """ 
    if str(type(f))=="<class 'function'>":
        #operator f je zadan formulom
        M=operator_kan(f)
    else:
        #operator f je zadan matricom u paru kanonskih baza
        M=f
    S=matrix(baza1).transpose()
    if S.det()==0: return "Error: prvi skup nije baza"
    if baza2==None:
        T=matrix(baza1).transpose()
    else:
        T=matrix(baza2).transpose()
    if T.det()==0: return "Error: drugi skup nije baza"
    return T^-1*M*S


def operator_baza(f,bazaD1,bazaD2,bazaK1=None,bazaK2=None):
    """ Ova funkcija je poopcenje funkcije operator_kanbaza. Ona daje
    matricni prikaz linearnog operatora f u paru baza (bazaD2,bazaK2)
    ako je on zadan u paru baza (bazaD1,bazaK1).
    Ako bazaK1 nije eksplicitno navedena, podrazumijeva se da je bazaK1=bazaD1.
    Ako bazaK2 nije eksplicitno navedena, podrazumijeva se da je bazaK2=bazaD2.
    
    Linearni operator f moze biti:
      * prethodno definiran kao funkcija u paru baza (bazaD1,bazaK1)
      * definiran na ulazu kao lambda funkcija u paru baza (bazaD1,bazaK1)
      * definiran svojim matricnim zapisom u paru baza (bazaD1,bazaK1)
    
    --- 1. PRIMJER ---
    
    Neka je [[2,0],[-1,0]] matrica operatora u paru baza {(1,1),(1,0)} i {(1,0),(0,1)}.
    Trazimo matricu tog operatora u bazi B={(1,0),(1,1)}.
    Dakle, u ovom slucaju je bazaD1={(1,1),(1,0)}, bazaD2=B, bazaK1={(1,0),(0,1)}, bazaK2=B.
    U ovom slucaju je bazaK2=bazaD2 pa ju nema potrebe dolje pisati na zadnjem mjestu.
    
    operator_baza(matrix([[2,0],[-1,0]]), [(1,1),(1,0)], [(1,0),(1,1)], [(1,0),(0,1)])
    operator_baza(matrix([[2,0],[-1,0]]), bazaD1=[(1,1),(1,0)], bazaD2=[(1,0),(1,1)], bazaK1=[(1,0),(0,1)])
    
    --- 2. PRIMJER ---
    
    Neka je [[1,0,2],[1,1,1],[1,0,1]] matrica linearnog operatora u bazi {(1,1,1),(1,2,1),(0,1,2)}.
    Trazimo matricu tog operatora u bazi {(1,0,1),(-1,0,1),(0,2,0)}.
    U ovom slucaju je bazaD1={(1,1,1),(1,2,1),(0,1,2)}, bazaD2={(1,0,1),(-1,0,1),(0,2,0)}, bazaK1=bazaD1, bazaK2=bazaD2.
    Kako je bazaK1=bazaD1, bazaK2=bazaD2, nema potrebe posebno pisati baze bazaK1, bazaK2 kod poziva funkcije.
    
    operator_baza(matrix([[1,0,2],[1,1,1],[1,0,1]]), [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)])
    
    Nije greska ako ih i posebno naglasimo.
    
    operator_baza(matrix([[1,0,2],[1,1,1],[1,0,1]]), [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)], [(1,1,1),(1,2,1),(0,1,2)], [(1,0,1),(-1,0,1),(0,2,0)])
    
    ili da nam bude preglednije
    
    operator_baza(matrix([[1,0,2],[1,1,1],[1,0,1]]), bazaD1=[(1,1,1),(1,2,1),(0,1,2)], bazaD2=[(1,0,1),(-1,0,1),(0,2,0)])
    operator_baza(matrix([[1,0,2],[1,1,1],[1,0,1]]), bazaD1=[(1,1,1),(1,2,1),(0,1,2)], bazaD2=[(1,0,1),(-1,0,1),(0,2,0)], bazaK1=[(1,1,1),(1,2,1),(0,1,2)], bazaK2=[(1,0,1),(-1,0,1),(0,2,0)])
    """     
    if str(type(f))=="<class 'function'>":
        #operator f je zadan formulom
        M=operator_kan(f)
    else:
        #operator f je zadan matricom u paru baza (bazaD1,bazaK1)
        M=f
    if bazaK1==None: bazaK1=bazaD1
    if bazaK2==None: bazaK2=bazaD2
    S=matrica_prijelaza(bazaD1,bazaD2)
    T=matrica_prijelaza(bazaK1,bazaK2)
    return T^-1*M*S


#### 5. poglavlje ####

def rucni_Horner(p,tocke):
    """Racuna vrijednost polinoma p u tocki x0 i na izlazu daje tablicu kakvu
    dobivamo rucnim izvodjenjem Hornerovog algoritma. 
    Varijabla p mora biti tipa 'SAGE polinom'. Koristite prstene QQ i SR za
    definiciju polinoma nad prstenom.
    
    Varijabla tocke moze biti samo jedan broj x0 ili pak lista brojeva ukoliko
    zelimo Hornera redom uzastopno primjenjivati na brojeve iz liste na
    prethodno dobiveni kvocijent.
    
    --- PRIMJERI ---
    
    P.<y>=QQ[]
    f1=4*y^3+2*y^2+y-1
    f2=8*y^4-76*y^3+270*y^2-425*y+250
    rucni_Horner(f1,5)
    rucni_Horner(f2,[5/2,5/2,5/2])
    R.<x>=SR[]
    g=(-1+I)*x^4+(1+I)*x^3+I*x+1
    f=x^6-3*x^5+12*x^4-21*x^3+47*x^2-36*x+60
    rucni_Horner(f,[2*I,-2*I,sqrt(3)*I,-sqrt(3)*I])
    rucni_Horner(g,I)
    """
    if str(type(tocke))!="<class 'list'>":
        tocke=[tocke]
    koef=list(reversed(p.coefficients()))
    tablica=[[" "]+koef]
    for redak in tocke:
        vrijednosti=[koef[0]]
        for k in range(1,len(koef)):
            vrijednosti.append(redak*vrijednosti[k-1]+koef[k])
        vrijednosti=list(map(pojednostavni,vrijednosti))
        tablica.append([redak]+vrijednosti)
        koef=vrijednosti[:-1]
    return table(tablica)


def razvoj_Horner(p,x0):
    """Daje razvoj polinoma p po potencijama od x-x0, tj. daje
    nam tablicu koju dobivamo rucnim rjesavanjem zadatka pomocu
    Hornerovog algoritma.
    Varijabla p mora biti tipa 'SAGE polinom'. Koristite prstene QQ i SR za
    definiciju polinoma nad prstenom.
    
    --- PRIMJER ---
    
    R.<x>=QQ[]
    f=4*x^3+2*x^2+x-1
    razvoj_Horner(f,-2)
    P.<y>=SR[]
    g=(-1+I)*y^4+(1+I)*y^3+I*y+1
    razvoj_Horner(g,2+sqrt(3)*I)
    """
    koef=list(reversed(p.coefficients()))
    d=len(koef)
    tablica=[[" "]+koef]
    for redak in range(d):
        vrijednosti=[koef[0]]
        for k in range(1,len(koef)):
            vrijednosti.append(x0*vrijednosti[k-1]+koef[k])
        vrijednosti=list(map(pojednostavni,vrijednosti))
        tablica.append([x0]+vrijednosti)
        koef=vrijednosti[:-1]
    return table(tablica)


def polinomi_gcd(f,g):
    """Daje svaki korak Euklidovog algoritma kod odredjivanja
    najvece zajednicke mjere polinoma f i g. U svakom koraku
    funkcija daje kvocijent i ostatak koji se dobiju dijeljenjem 
    odgovarajucih polinoma.
    Varijable f i g moraju biti tipa 'SAGE polinom'. Koristite prstene 
    QQ i SR za definiciju polinoma nad prstenom.
    
    --- PRIMJER ---
    
    R.<x>=QQ[]
    polinomi_gcd(x^5+x+1, x^4+2*x^3+1)
    """
    tablica=[['korak','------------------- kvocijent -------------------','-------------------  ostatak  -------------------']]
    korak=1
    q,r=f,g
    while r!=0:
        tablica.append([str(korak)+'.',q//r,q%r])
        q,r=r,q%r
        korak += 1
    return table(tablica)


def kompleksni_korijen(z,n,sredi=False):
    """Daje egzaktno sve n-te korijene kompleksnog broja z.
    
    sredi=False -> nece se sredjivati kosinusi i sinusi lijepih kutova
    sredi=True  -> sredit ce se kosinusi i sinusi lijepih kutova

    --- PRIMJERI ---
    
    kompleksni_korijen(9,5)
    kompleksni_korijen(9,5,sredi=True)
    
    kompleksni_korijen(-2-2*I,7)
    
    kompleksni_korijen(-1/2+sqrt(3)/2*I,3)
    kompleksni_korijen(-1/2+sqrt(3)/2*I,3,sredi=True)
    """
    a,b=z.real(),z.imag()
    r=simplify(abs(z))
    if a==0 and b==0: return 0
    if a==0 and b>0: fi=pi/2
    if a==0 and b<0: fi=3*pi/2
    if a>0 and b>=0: fi=simplify(atan(b/a))
    if a<0 and b>=0: fi=simplify(atan(b/a))+pi
    if a<0 and b<0: fi=simplify(atan(b/a))+pi
    if a>0 and b<0: fi=simplify(atan(b/a))+2*pi
    kutovi=[(fi+2*k*pi)/n for k in range(n)]
    if sredi:
        korijeni=list(map(lambda t: factor(r^(1/n)*(cos(t,hold=True)+I*sin(t,hold=True))),kutovi))
    else:
        korijeni=list(map(lambda t: r^(1/n)*(cos(t,hold=True)+I*sin(t,hold=True)),kutovi))
    return korijeni


def Ferrari_rezolventa(p):
    """Daje Ferrarijevu rezolventu algebarske jednadzbe 4. stupnja.
    Varijabla p mora biti tipa 'SAGE polinom', a na izlazu se vraca 
    Ferrarijeva rezolventa u varijabli y. Koristite prstene 
    QQ i SR za definiciju polinoma nad prstenom.
    
    --- PRIMJER ---
    
    R.<x>=QQ[]
    Ferrari_rezolventa(x^4+3*x^2+2*x+3)
    """
    koef=p.coefficients(sparse=False)
    if len(koef)!=5: return "Error: polinom mora biti 4. stupnja"
    a,b,c,d=koef[3]/koef[4],koef[2]/koef[4],koef[1]/koef[4],koef[0]/koef[4]
    y=var('y')
    izraz=expand((a*y-c)^2-4*(a^2/4-b+2*y)*(y^2-d))
    return izraz==0


def Cardan(p,q):
    """Daje sva rjesenja jednadzbe x^3+px+q=0 koristeci Cardanovu formulu.
    
    --- PRIMJER ---
    
    Cardan(-3,1) #rjesenja jednadzbe x^3-3x+1=0
    Cardan(-12,-8) #rjesenja jednadzbe x^3-12x-8=0
    """
    a1=(q/2)^2+(p/3)^3
    if a1>=0:
        a2=-q/2+sqrt(a1)
    else:
        a2=-q/2+sqrt(-a1)*I
    if a2.imag()==0:
        u=a2^(1/3)
    else:
        u=kompleksni_korijen(a2,3,sredi=True)[0]
    v=expand(-p/(3*u))
    if v.imag()!=0:
        v=(-p/(3*u)).real().trig_simplify()+(-p/(3*u)).imag().trig_simplify()*I
    e1=-1/2-sqrt(3)/2*I
    e2=-1/2+sqrt(3)/2*I
    x1=expand(u+v)
    x2=expand(u*e1+v*e2)
    x3=expand(u*e2+v*e1)
    return [x1,x2,x3]
