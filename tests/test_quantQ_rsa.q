/ Unit tests for lib/quantQ_rsa.q -- number theory and textbook RSA (.quantQ.rsa)
/ Reference values from Python (pow(e,-1,n), pow(b,e,n)).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`rsa;

/ primes
.qt.eq["isPrime on primes";{.quantQ.rsa.isPrime each 2 3 5 7 31 97 127};7#1b];
.qt.eq["isPrime on composites";{.quantQ.rsa.isPrime each 4 9 15 21 91 121};6#0b];
.qt.eq["primeUp finds the next prime";{.quantQ.rsa.primeUp 14};17];
.qt.eq["primeUp of a prime is itself";{.quantQ.rsa.primeUp 17};17];
.qt.eq["primeUp used by genKeys";{.quantQ.rsa.primeUp each 123 345};127 347];

/ gcd / lcm / coprimality
.qt.eq["gcd";{.quantQ.rsa.gcd[14;21]};7];
.qt.eq["gcd is symmetric";{.quantQ.rsa.gcd[21;14]};7];
.qt.eq["gcd of coprimes";{.quantQ.rsa.gcd[17;5]};1];
.qt.true["isCoprime";{.quantQ.rsa.isCoprime[14;15]}];
.qt.true["not coprime";{not .quantQ.rsa.isCoprime[14;28]}];
.qt.approx["lcm";{.quantQ.rsa.lcm[6;14]};42;0f];
.qt.approx["lcm of zeros is 0";{.quantQ.rsa.lcm[0;0]};0;0f];

/ Carmichael function lambda(n)
.qt.eq["CarmichaelFunc(14) = 6";{.quantQ.rsa.CarmichaelFunc 14};6];
.qt.eq["CarmichaelFunc(15) = 4";{.quantQ.rsa.CarmichaelFunc 15};4];
.qt.eq["CarmichaelFunc(8) = 2";{.quantQ.rsa.CarmichaelFunc 8};2];

/ modular arithmetic
.qt.eq["inverseMod: 17^-1 mod 21798";{.quantQ.rsa.inverseMod[17;21798]`t};5129];
.qt.eq["inverseMod: status 1 when the inverse exists";{.quantQ.rsa.inverseMod[17;21798]`status};1];
.qt.eq["inverseMod: status 0 when it does not";{.quantQ.rsa.inverseMod[2;4]`status};0];
.qt.eq["modExp: 4^13 mod 497";{.quantQ.rsa.modExp[4;13;497]};445];
.qt.eq["modExp: exponent 0";{.quantQ.rsa.modExp[4;0;497]};1];

/ keys, encryption and decryption
rsaKeys:.qt.fixture["rsaKeys";{.quantQ.rsa.genKeys ()!()}];
.qt.eq["genKeys: default primes";{rsaKeys`p`q};127 347];
.qt.eq["genKeys: modulus";{rsaKeys`modulus};44069];
.qt.eq["genKeys: lambda = lcm(p-1,q-1)";{rsaKeys`lam};21798];
.qt.eq["genKeys: public exponent";{rsaKeys`e};17];
.qt.eq["genKeys: private exponent";{rsaKeys`d};5129];
.qt.eq["genKeys: public key";{rsaKeys`publicKey};`modulus`e!44069 17];
.qt.eq["encrypt";{.quantQ.rsa.encrypt[rsaKeys`publicKey;123]`cipher};3111];
.qt.eq["encrypt refuses messages >= modulus";{.quantQ.rsa.encrypt[rsaKeys`publicKey;50000]`cipher};0];
.qt.eq["decrypt(encrypt(m)) = m";{.quantQ.rsa.decrypt[.quantQ.rsa.encrypt[rsaKeys`publicKey;123];rsaKeys]};123];
.qt.eq["round trip for several messages";{{.quantQ.rsa.decrypt[.quantQ.rsa.encrypt[rsaKeys`publicKey;x];rsaKeys]} each 2 99 1000 44068};2 99 1000 44068];

.qt.finish[];
