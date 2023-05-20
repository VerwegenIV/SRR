\ Model Obj
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  _var2_ + _var3_ + _var5_ + _var6_ + _var7_ + 3 _var8_ + _var9_ + 2 _var11_
   + 2 _var12_ + 2 _var13_ + _var14_ + 2 _var17_ + 3 _var18_ + 2 _var19_
   + _var20_ + _var21_ + 3 _var22_ + _var24_ + 2 _var25_ + 2 _var26_
   + 3 _var27_ + _var28_ + _var29_ + _var30_ + 2 _var31_ + 2 _var32_
   + 2 _var33_ + _var35_ + 2 _var36_ + _var37_ + _var38_ + 3 _var39_
   + 2 _var40_ + _var41_ + 2 _var42_ + _var43_ + 3 _var44_ + 2 _var46_
   + 2 _var47_ + 2 _var48_ + _var49_ + 2 _var50_ + 2 _var51_ + 2 _var53_
Subject To
 round[0]: _var3_ + _var6_ + _var9_ + _var14_ + _var19_ + _var24_ + _var29_
   + _var39_ + _var44_ + _var0_ + _var34_ = 1
 round[1]: _var7_ + _var20_ + _var25_ + _var30_ + _var35_ + _var40_
   + _var49_ + _var1_ + _var4_ + _var10_ + _var15_ + _var45_ = 1
 round[2]: _var2_ + _var5_ + _var8_ + _var11_ + _var21_ + _var26_ + _var31_
   + _var36_ + _var41_ + _var46_ + _var50_ + _var16_ = 1
 round[3]: _var12_ + _var17_ + _var22_ + _var27_ + _var32_ + _var37_
   + _var42_ + _var47_ + _var51_ + _var53_ = 1
 round[4]: _var13_ + _var18_ + _var28_ + _var33_ + _var38_ + _var43_
   + _var48_ + _var23_ + _var52_ = 1
 match[0: _var2_ + _var6_ + _var7_ + _var8_ + _var24_ + _var25_ + _var26_
   + _var27_ + _var28_ + _var0_ + _var1_ = 1
 match[1: _var3_ + _var5_ + _var35_ + _var36_ + _var37_ + _var38_ + _var34_
   + _var4_ + _var52_ = 1
 match[2: _var14_ + _var17_ + _var18_ + _var39_ + _var40_ + _var41_
   + _var42_ + _var43_ + _var15_ + _var16_ = 1
 match[3: _var9_ + _var11_ + _var12_ + _var13_ + _var29_ + _var30_
   + _var31_ + _var32_ + _var33_ + _var49_ + _var50_ + _var51_ + _var10_
   + _var53_ = 1
 match[4: _var19_ + _var20_ + _var21_ + _var22_ + _var44_ + _var46_
   + _var47_ + _var48_ + _var45_ + _var23_ = 1
 match[5: _var14_ + _var17_ + _var18_ + _var29_ + _var30_ + _var31_
   + _var32_ + _var33_ + _var15_ + _var16_ = 1
 match[6: _var3_ + _var5_ + _var44_ + _var46_ + _var47_ + _var48_ + _var49_
   + _var50_ + _var51_ + _var4_ + _var45_ = 1
 match[7: _var19_ + _var20_ + _var21_ + _var22_ + _var39_ + _var40_
   + _var41_ + _var42_ + _var43_ + _var23_ + _var52_ = 1
 match[8: _var9_ + _var11_ + _var12_ + _var13_ + _var35_ + _var36_
   + _var37_ + _var38_ + _var34_ + _var10_ + _var53_ = 1
 match[9: _var2_ + _var9_ + _var11_ + _var12_ + _var13_ + _var19_ + _var20_
   + _var21_ + _var22_ + _var0_ + _var1_ + _var10_ + _var23_ + _var53_ = 1
 match[10: _var6_ + _var7_ + _var8_ + _var44_ + _var46_ + _var47_ + _var48_
   + _var45_ = 1
 match[11: _var24_ + _var25_ + _var26_ + _var27_ + _var28_ + _var39_
   + _var40_ + _var41_ + _var42_ + _var43_ + _var49_ + _var50_ + _var51_
   = 1
 match[12: _var24_ + _var25_ + _var26_ + _var27_ + _var28_ + _var35_
   + _var36_ + _var37_ + _var38_ + _var34_ = 1
 match[13: _var6_ + _var7_ + _var8_ + _var29_ + _var30_ + _var31_ + _var32_
   + _var33_ + _var52_ = 1
 match[14: _var2_ + _var3_ + _var5_ + _var14_ + _var17_ + _var18_ + _var0_
   + _var1_ + _var4_ + _var15_ + _var16_ = 1
 cut0: - _var9_ - _var11_ - _var13_ + _var17_ + _var22_ + _var27_ + _var42_
   + _var47_ - _var10_ >= 0
 cut1: _var7_ - _var19_ - _var21_ - _var22_ + _var25_ + _var30_ + _var35_
   + _var40_ + _var49_ + _var4_ + _var15_ - _var23_ >= 0
 cut2: _var12_ + _var17_ + _var22_ - _var24_ - _var25_ - _var26_ - _var28_
   + _var32_ + _var47_ + _var53_ >= 0
 cut3: _var7_ + _var20_ + _var25_ + _var35_ + _var40_ - _var50_ - _var51_
   + _var1_ + _var15_ >= 0
Bounds
 _var2_ <= 1
 _var3_ <= 1
 _var5_ <= 1
 _var6_ <= 1
 _var7_ <= 1
 _var8_ <= 1
 _var9_ <= 1
 _var11_ <= 1
 _var12_ <= 1
 _var13_ <= 1
 _var14_ <= 1
 _var17_ <= 1
 _var18_ <= 1
 _var19_ <= 1
 _var20_ <= 1
 _var21_ <= 1
 _var22_ <= 1
 _var24_ <= 1
 _var25_ <= 1
 _var26_ <= 1
 _var27_ <= 1
 _var28_ <= 1
 _var29_ <= 1
 _var30_ <= 1
 _var31_ <= 1
 _var32_ <= 1
 _var33_ <= 1
 _var35_ <= 1
 _var36_ <= 1
 _var37_ <= 1
 _var38_ <= 1
 _var39_ <= 1
 _var40_ <= 1
 _var41_ <= 1
 _var42_ <= 1
 _var43_ <= 1
 _var44_ <= 1
 _var46_ <= 1
 _var47_ <= 1
 _var48_ <= 1
 _var49_ <= 1
 _var50_ <= 1
 _var51_ <= 1
 _var0_ <= 1
 _var34_ <= 1
 _var1_ <= 1
 _var4_ <= 1
 _var10_ <= 1
 _var15_ <= 1
 _var45_ <= 1
 _var16_ <= 1
 _var23_ <= 1
 _var52_ <= 1
 _var53_ <= 1
End
