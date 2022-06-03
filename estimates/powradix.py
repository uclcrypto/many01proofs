from gmpy2 import mpz, xmpz, mpz_urandomb, powmod, random_state
from timeit import timeit

# Number of exponentiations to be computed in a single basis
n = 1000
# Size of the exponent
e_size = 256

# 3072 bits p
# p = mpz(1625342247057486122165665607533780094243597738493809630880578236306193245103134066809379590865889503357217507435006216465066227850714707580196090135702171306627525077342334624315728843216085190499113122280899917450884661126702195194345733278270642438105592610410552260717002903549487058012322754637448658162339576786827461787115999291079065111490635596424716216903658564885813357260817414492342592852345984145354034039825266919075610537869394654501240956613594199942228090535091043991889817748340728690806845964043784420056083818626719317961681516642950882927723197601986900120741850785026878262645352821297276511462388355980067252357682607298913819383011033767798795353308117011977881858244274869398520064145364602126321305994802172725168359767468742681333010031895125325095132032571245228824338971152948569286032532599959795301053964807472092382975399104660938652873450290742589627747591394267203993896260119571790772029339)
# g = mpz(2)

# ElectionGuard 4096 bits p, and g
p = mpz("0x FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 93C467E3 7DB0C7A4 D1BE3F81 0152CB56 A1CECC3A F65CC019 0C03DF34 709AFFBD 8E4B59FA 03A9F0EE D0649CCB 621057D1 1056AE91 32135A08 E43B4673 D74BAFEA 58DEB878 CC86D733 DBE7BF38 154B36CF 8A96D156 7899AAAE 0C09D4C8 B6B7B86F D2A1EA1D E62FF864 3EC7C271 82797722 5E6AC2F0 BD61C746 961542A3 CE3BEA5D B54FE70E 63E6D09F 8FC28658 E80567A4 7CFDE60E E741E5D8 5A7BD469 31CED822 03655949 64B83989 6FCAABCC C9B31959 C083F22A D3EE591C 32FAB2C7 448F2A05 7DB2DB49 EE52E018 2741E538 65F004CC 8E704B7C 5C40BF30 4C4D8C4F 13EDF604 7C555302 D2238D8C E11DF242 4F1B66C2 C5D238D0 744DB679 AF289048 7031F9C0 AEA1C4BB 6FE9554E E528FDF1 B05E5B25 6223B2F0 9215F371 9F9C7CCC 69DDF172 D0D62342 17FCC003 7F18B93E F5389130 B7A661E5 C26E5421 4068BBCA FEA32A67 818BD307 5AD1F5C7 E9CC3D17 37FB2817 1BAF84DB B6612B78 81C1A48E 439CD03A 92BF5222 5A2B38E6 542E9F72 2BCE15A3 81B5753E A8427633 81CCAE83 512B3051 1B32E5E8 D8036214 9AD030AA BA5F3A57 98BB22AA 7EC1B6D0 F17903F4 E234EA60 34AA8597 3F79A93F FB82A75C 47C03D43 D2F9CA02 D03199BA CEDDD453 34DBC6B5 FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF")
g = mpz("0x 037DE384 F98F6E03 8D2A3141 825B33D5 D45EC4CC 64CFD15E 750D6798 F5196CF2 A142CDF3 3F6EF853 840EC7D4 EC804794 CFB0CFB6 5363B256 6387B98E E0E3DEF1 B706FA55 D5038FFB 4A62DCBB 93B1DDD8 D3B308DA 86D1C3A5 25EF356F E5BB5931 4E656334 80B396E1 DD4B795F 78DE07D8 6B0E2A05 BE6AF78F D7F736FC BA6C032E 26E050AF 50A03C65 FA7B6C87 F4554CB5 7F3DABCB AD8EB9D8 FDEBEEF5 8570669A CC3EDA17 DBFC47B8 B3C39AA0 8B829B28 872E62B5 D1B13A98 F09D40AC 20C2AB74 A6750E7C 8750B514 1E221C41 F55BBA31 D8E41422 B64D2CBA 7AAA0E9F D8785702 F6932825 BF45DE83 86D24900 742062C1 322B37C5 0AF18215 8090C35D A9355E6C F7F72DA3 9A2284FD FB1918B2 A2A30E69 501FA234 2B728263 DF23F1DB 8355BDE1 EB276FB3 685F3716 72CEB313 FDAB069C C9B11AB6 C59BCE62 BAAD96AA C96B0DBE 0C7E71FC B2255254 5A5D1CED EEE01E4B C0CDBDB7 6B6AD45F 09AF5E71 114A005F 93AD97B8 FE09274E 76C94B20 08926B38 CAEC94C9 5E96D628 F6BC8066 2BA06207 801328B2 C6A60526 BF7CD02D 9661385A C3B1CBDB 50F759D0 E9F61C11 A07BF421 8F299BCB 29005200 76EBD2D9 5A3DEE96 D4809EF3 4ABEB83F DBA8A12C 5CA82757 288A89C9 31CF564F 00E8A317 AE1E1D82 8E61369B A0DDBADB 10C136F8 691101AD 82DC5477 5AB83538 40D99921 97D80A6E 94B38AC4 17CDDF40 B0C73ABF 03E8E0AA")


# Picking a list of n exponents
seed = random_state()
e_list = [xmpz(mpz_urandomb(seed, e_size)) for i in range(n)]


# Basic square and multiply, precomputing squares, and taking advantage of iterations on xmpz
# It is equivalent to PowRadix for k=1 but runs slightly faster
class PowRadix2:
    def __init__(self, basis):
        squares = []
        gs = basis
        for i in range(e_size):
            squares.append(gs)
            gs = gs * gs % p
        self.squares = squares

    def pow(self, e):
        y = mpz(1)
        for i in e.iter_set():
            y = y * self.squares[i] % p
        return y


# Radix method
class PowRadix:
    def __init__(self, basis, k=1, n=None):
        # if n is given, then looking for the best k
        if n:
            k = 1
            while (0.69*k-1) * 2**k < n:  # Equality happens for optimal k
                k += 1
            k -= 1  # limiting amount of precomputation
        self.table_length = -(-e_size // k)  # Double negative to take the ceiling
        self.k = k
        table = []
        row_basis = basis
        running_basis = row_basis
        for _ in range(self.table_length):
            row = [1]
            for j in range(1, 2**k):
                row.append(running_basis)
                running_basis = running_basis * row_basis % p
            table.append(row)
            row_basis = running_basis
        self.table = table

    def pow(self, e):
        y = mpz(1)
        for i in range(self.table_length):
            e_slice = e[i * self.k: (i+1) * self.k]
            y = y * self.table[i][e_slice] % p
        return y

    def alt_pow(self, e):
        # Trying to see if this runs faster, but it does not
        y = mpz(1)
        slice_start = 0
        for row in self.table:
            slice_end = slice_start + self.k
            e_slice = e[slice_start: slice_end]
            slice_start = slice_end
            y = y * row[e_slice] % p
        return y


# Testing and initializing before benchmarks
yp = [powmod(g, e, p) for e in e_list]

g_radix_2 = PowRadix2(g)
y_radix_2 = [g_radix_2.pow(e) for e in e_list]
assert(yp == y_radix_2)

g_radix = PowRadix(g, n=n)
y_radix = [g_radix.pow(e) for e in e_list]
assert(yp == y_radix)

# Benchmarks
iterations = 10
t1 = timeit('[powmod(g, e, p) for e in e_list]', globals=globals(), number=iterations) * 1000 / iterations
print("Computing", n, "exponentiations with gmpy2: ", t1, "ms.")

t2 = timeit('PowRadix(g,n=n)', globals=globals(), number=iterations) * 1000 / iterations
print("Precomputation for radix method with k =", g_radix.k, ": ", t2, "ms.")

t3 = timeit('[g_radix.pow(e) for e in e_list]', globals=globals(), number=iterations) * 1000 / iterations
print("Computing", n, "exponentiations with radix method: ", t3, "ms.")

print("Speedup with radix method: ", t1/(t2+t3))

t4 = timeit('PowRadix2(g)', globals=globals(), number=iterations) * 100
print("Precomputation for radix 2 method: ", t4, "ms.")

t5 = timeit('[g_radix_2.pow(e) for e in e_list]', globals=globals(), number=iterations) * 1000 / iterations
print("Computing", n, "exponentiations with radix 2 method: ", t5, "ms.")

print("Speedup with radix 2 method: ", t1/(t4+t5))
