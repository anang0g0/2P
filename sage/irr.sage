BIN.<X> = GF(16,"a")[]
while 1:
      poly = BIN.random_element(2);
      if poly.is_irreducible():
      	 break;
print(poly)
 
