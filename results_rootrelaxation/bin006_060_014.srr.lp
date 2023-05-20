\ SCIP STATISTICS
\   Problem name     : t_instances/bin006_060_014.srr
\   Variables        : 54 (0 binary, 0 integer, 0 implicit integer, 54 continuous)
\   Constraints      : 20
Minimize
 Obj: +2 _var0_ +2 _var1_ +2 _var2_ +3 _var3_ +3 _var4_ +2 _var5_ +2 _var8_ +2 _var9_ +3 _var10_ +2 _var11_
      +3 _var12_ +3 _var13_ +2 _var14_ +1 _var15_ +2 _var16_ +3 _var17_ +1 _var18_ +2 _var19_ +2 _var20_
      +1 _var21_ +3 _var22_ +2 _var23_ +1 _var24_ +1 _var25_ +2 _var26_ +1 _var27_ +2 _var28_ +1 _var30_
      +3 _var31_ +2 _var32_ +1 _var33_ +2 _var34_ +3 _var35_ +2 _var36_ +3 _var37_ +2 _var38_ +3 _var39_
      +2 _var42_ +1 _var43_ +1 _var44_ +2 _var45_ +2 _var46_ +3 _var47_ +1 _var48_ +1 _var50_ +1 _var51_
      +2 _var52_ +1 _var53_
Subject to
 round[0]: +1 _var0_ +1 _var3_ +1 _var6_ +1 _var9_ +1 _var14_ +1 _var19_ +1 _var24_ +1 _var29_ +1 _var34_
  +1 _var39_ +1 _var44_ +1 _var49_ = +1
 round[1]: +1 _var1_ +1 _var4_ +1 _var7_ +1 _var10_ +1 _var15_ +1 _var20_ +1 _var25_ +1 _var30_ +1 _var35_
  +1 _var40_ +1 _var45_ +1 _var50_ = +1
 round[2]: +1 _var2_ +1 _var5_ +1 _var8_ +1 _var11_ +1 _var16_ +1 _var21_ +1 _var26_ +1 _var31_ +1 _var36_
  +1 _var41_ +1 _var46_ +1 _var51_ = +1
 round[3]: +1 _var12_ +1 _var17_ +1 _var22_ +1 _var27_ +1 _var32_ +1 _var37_ +1 _var42_ +1 _var47_ +1 _var52_
  = +1
 round[4]: +1 _var13_ +1 _var18_ +1 _var23_ +1 _var28_ +1 _var33_ +1 _var38_ +1 _var43_ +1 _var48_ +1 _var53_
  = +1
 match[0:0,1]: +1 _var0_ +1 _var1_ +1 _var2_ +1 _var6_ +1 _var7_ +1 _var8_ +1 _var24_ +1 _var25_ +1 _var26_
  +1 _var27_ +1 _var28_ +1 _var52_ = +1
 match[1:0,2]: +1 _var3_ +1 _var4_ +1 _var5_ +1 _var34_ +1 _var35_ +1 _var36_ +1 _var37_ +1 _var38_ = +1
 match[2:0,3]: +1 _var14_ +1 _var15_ +1 _var16_ +1 _var17_ +1 _var18_ +1 _var39_ +1 _var40_ +1 _var41_
  +1 _var42_ +1 _var43_ +1 _var50_ +1 _var51_ +1 _var53_ = +1
 match[3:0,4]: +1 _var9_ +1 _var10_ +1 _var11_ +1 _var12_ +1 _var13_ +1 _var29_ +1 _var30_ +1 _var31_
  +1 _var32_ +1 _var33_ = +1
 match[4:0,5]: +1 _var19_ +1 _var20_ +1 _var21_ +1 _var22_ +1 _var23_ +1 _var44_ +1 _var45_ +1 _var46_
  +1 _var47_ +1 _var48_ +1 _var49_ = +1
 match[5:1,2]: +1 _var14_ +1 _var15_ +1 _var16_ +1 _var17_ +1 _var18_ +1 _var29_ +1 _var30_ +1 _var31_
  +1 _var32_ +1 _var33_ +1 _var49_ = +1
 match[6:1,3]: +1 _var3_ +1 _var4_ +1 _var5_ +1 _var44_ +1 _var45_ +1 _var46_ +1 _var47_ +1 _var48_ = +1
 match[7:1,4]: +1 _var19_ +1 _var20_ +1 _var21_ +1 _var22_ +1 _var23_ +1 _var39_ +1 _var40_ +1 _var41_
  +1 _var42_ +1 _var43_ = +1
 match[8:1,5]: +1 _var9_ +1 _var10_ +1 _var11_ +1 _var12_ +1 _var13_ +1 _var34_ +1 _var35_ +1 _var36_
  +1 _var37_ +1 _var38_ +1 _var50_ +1 _var51_ +1 _var53_ = +1
 match[9:2,3]: +1 _var0_ +1 _var1_ +1 _var2_ +1 _var9_ +1 _var10_ +1 _var11_ +1 _var12_ +1 _var13_ +1 _var19_
  +1 _var20_ +1 _var21_ +1 _var22_ +1 _var23_ +1 _var52_ = +1
 match[10:2,4]: +1 _var6_ +1 _var7_ +1 _var8_ +1 _var44_ +1 _var45_ +1 _var46_ +1 _var47_ +1 _var48_ +1 _var50_
  +1 _var51_ +1 _var53_ = +1
 match[11:2,5]: +1 _var24_ +1 _var25_ +1 _var26_ +1 _var27_ +1 _var28_ +1 _var39_ +1 _var40_ +1 _var41_
  +1 _var42_ +1 _var43_ = +1
 match[12:3,4]: +1 _var24_ +1 _var25_ +1 _var26_ +1 _var27_ +1 _var28_ +1 _var34_ +1 _var35_ +1 _var36_
  +1 _var37_ +1 _var38_ +1 _var49_ = +1
 match[13:3,5]: +1 _var6_ +1 _var7_ +1 _var8_ +1 _var29_ +1 _var30_ +1 _var31_ +1 _var32_ +1 _var33_ = +1
 match[14:4,5]: +1 _var0_ +1 _var1_ +1 _var2_ +1 _var3_ +1 _var4_ +1 _var5_ +1 _var14_ +1 _var15_ +1 _var16_
  +1 _var17_ +1 _var18_ +1 _var52_ = +1
Bounds
 0 <= _var0_ <= 1
 0 <= _var1_ <= 1
 0 <= _var2_ <= 1
 0 <= _var3_ <= 1
 0 <= _var4_ <= 1
 0 <= _var5_ <= 1
 0 <= _var6_ <= 1
 0 <= _var7_ <= 1
 0 <= _var8_ <= 1
 0 <= _var9_ <= 1
 0 <= _var10_ <= 1
 0 <= _var11_ <= 1
 0 <= _var12_ <= 1
 0 <= _var13_ <= 1
 0 <= _var14_ <= 1
 0 <= _var15_ <= 1
 0 <= _var16_ <= 1
 0 <= _var17_ <= 1
 0 <= _var18_ <= 1
 0 <= _var19_ <= 1
 0 <= _var20_ <= 1
 0 <= _var21_ <= 1
 0 <= _var22_ <= 1
 0 <= _var23_ <= 1
 0 <= _var24_ <= 1
 0 <= _var25_ <= 1
 0 <= _var26_ <= 1
 0 <= _var27_ <= 1
 0 <= _var28_ <= 1
 0 <= _var29_ <= 1
 0 <= _var30_ <= 1
 0 <= _var31_ <= 1
 0 <= _var32_ <= 1
 0 <= _var33_ <= 1
 0 <= _var34_ <= 1
 0 <= _var35_ <= 1
 0 <= _var36_ <= 1
 0 <= _var37_ <= 1
 0 <= _var38_ <= 1
 0 <= _var39_ <= 1
 0 <= _var40_ <= 1
 0 <= _var41_ <= 1
 0 <= _var42_ <= 1
 0 <= _var43_ <= 1
 0 <= _var44_ <= 1
 0 <= _var45_ <= 1
 0 <= _var46_ <= 1
 0 <= _var47_ <= 1
 0 <= _var48_ <= 1
 0 <= _var49_ <= 1
 0 <= _var50_ <= 1
 0 <= _var51_ <= 1
 0 <= _var52_ <= 1
 0 <= _var53_ <= 1
End
