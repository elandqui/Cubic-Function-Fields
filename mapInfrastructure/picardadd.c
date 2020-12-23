#define P 7 /* The characteristic. */

/* Globals */
int res1;
int u1[3] = {2,0,2};
int u2[3] = {0,3,2};
int v1[3] = {0,1,0};
int v2[3] = {1,3,1};
int f[5] = {1, 2,1,0,1}; /* y^3 = f[0] + f[1]x + ... + f[4]x^4. */
int inv0, inv1, inv2; /* z1 = inv0 x^2 + inv1 x + inv2 */
int t62, t63, t64, t65, t66, t67; /* s = t62 x + t63, t = t64 x^3 + ... + t67 */
int vsub[3];
int c1, c2, c3;

int inverse(a) {
  if(a==2) return 4;
  else if(a==3) return 5;
  else if(a==4) return 2;
  else if(a==5) return 3;
  else return a;
}

void resultant(){
  int t1 = u2[1]*vsub[2] %P;
  int t2 = u2[2]*vsub[2] %P;
  int t3 = u2[0]*vsub[2] %P;
  int t4 = u2[2]*vsub[0] %P;
  int t5 = u2[1]*vsub[1] %P;
  int t6 = vsub[2]*(t1+vsub[0]) %P;
  int t7 = vsub[1]*(vsub[1]-t2) %P;
  int t8 = (t4-t3-t5)*(vsub[1]+t2) %P;
  int t9 = vsub[2]*(t4-t3-t5) %P;
  int t10 = vsub[1]*(vsub[0]-t1) %P;
  inv0 = t6+t7 %P;
  int t11 = inv0*u2[2] %P;
  int t12 = u2[0]*vsub[1] %P;
  /* int t13 = inv0*t12 %P; */
  /* int t14 = t3*(t9-t10) %P; */
  /* int s1 = (v[0]-t1)*(v[0]-t1)%P; */
  inv2 = t8+ (vsub[0]-t1)*(vsub[0]-t1)%P;
  /* int t15 = inv2*v[0] %P;*/
  inv1 = t11 + t9 - t10 %P;
  res1  = (inv2*vsub[0] %P) - (inv0*t12 %P) - (t3*(t9-t10) %P) %P;
}

/* Computes the cubic w = y^2 + sy + t. */
void cubic( ){
  int usub[3];
  usub[0] = u1[0]-u2[0];
  usub[1] = u1[1]-u2[1];
  usub[2] = u1[2]-u2[2];

  int t16 = usub[2]*inv0 %P;
  int t17 = usub[1]*inv1 %P;
  int t18 = usub[0]*inv2 %P;
  int t19 = (usub[1] + usub[2])*(inv0+inv1)%P;
  int t20 = (usub[0] + usub[2])*(inv0+inv2)%P;
  int t21 = (usub[0] + usub[1])*(inv2+inv1)%P;
  int t22 = u2[2]*t16 %P;
  int t23 = u2[1]*t16 %P;
  int t24 = u2[2]*(t22+t16+t17-t19) %P;
  int t25 = (u2[1] + u2[0] %P)*(t19-t22-t17%P)%P;
  int t26 = u2[0]*(t22+t16+t17-t19%P)%P;
  int r0 = t24 + t20 + t17 - t23 -t16 - t18 %P;
  int r1 = t21 + t23 - t17 - t18 - t25 - t26 %P;
  int r2 = t18 + t26 %P;
  int s2 = v1[2]*v1[2] %P;
  int t27 = r0*res1 %P;
  int t28 = r0*s2 %P;
  int t29 = r0*t28 %P;
  int t30 = t28*res1 %P;
  int t31 = -res1*(v1[2] + v2[2])%P;
  int t32 = r1*s2%P;
  int t33 = u2[2]*t28%P;
  int g1 = t31 + t33 -t32%P;
  int t34 = res1*g1 %P;
  int t35 = -t27*(v1[1]+v2[1])%P;
  int t36 = -t27*(v1[0] + v2[0])%P;
  int t37 = r1*g1 %P;
  int t38 = r2*t28%P;
  int t39 = r2*g1 %P;
  int t40 = u2[1]*t29 %P;
  int t41 = u2[0]*t29%P;
  int l1 = t35+t40-t37-t38%P;
  int m1  = t36+t41 - t39%P;
  int t42 = -t27*v1[2]%P;
  int t43 = -t27*v1[1]%P;
  int t44 = -t27*v1[0]%P;
  int t45 = (v1[2]+v1[1])*(t42+t43-l1)%P;
  int t46 = v1[1]*(t43-l1)%P;
  int t47 = (v1[2]+v1[0])*(t42+t44-m1)%P;
  int t48 = v1[0]*(t44-m1)%P;
  int t49 = (v1[1]+v1[0])*(t43+t44-l1-m1)%P;
  int t50 = t30*(u1[2]+u1[1])%P;
  int t51 = u1[1]*t30%P;
  int t52 = t34*(u1[2]+u1[0])%P;
  int t53 = u1[0]*t34 %P;
  int t54 = (u1[1]+u1[0])*(t30+t34)%P;
  int B0 = t34+t50+t45+t30 - t51-t46%P;
  int B1 = t52+t30+t51+t47+t46-t53-t48%P;
  int B2 = t54+t49-t51-t53-t46-t48%P;
  int B3 = t53 + t48 %P;
  int t55 = B0*t27%P;
  int i1 = inverse(t55);
  int t56 = i1*B0 %P;
  int t57 = i1*t27 %P;
  int t58 = t57*t27%P;
  int t59 = t57*B1%P;
  int t60 = t57*B2%P;
  int t61 = t57*B3%P;
  t62 = t56*l1%P;
  t63 = t56*m1%P;
  t64 = t56*B0%P;
  t65 = t56*B1%P;
  t66 = t56*B2%P;
  t67 = t56*B3%P;

  /* Compute res(w,C,y) */
  int s3 = t59*t59%P;
  int t68 = t59*(6*t60+s3)%P;
  int s4 = t62*t62%P;
  int s5 = (t62+t63)*(t62+t63)%P;
  int s6 = t63*t63%P;
  int t69 = t62*t64%P;
  int t70 = t62*(s4-3*t65)%P;
  int t71 = t63*t64%P;
  int t72 = -3*f[3]*t69%P;
  int t73 = t62*(s5-3*t66-s4-s6)%P;
  int t74 = t63*(s4-3*t65)%P;
  int t75 = f[3]*t70%P;
  int t76 = -3*f[2]*t69%P;
  int t77 = -3*f[3]*t71%P;
  int s7 = t58*t58%P;
  int t78 = t58*s7%P;
  int t79 = t78*(1-3*t69)%P;
  int t80 = t78*(t70+t72+2*f[3]-3*t71)%P;
  int t81 = t78*(t73+t74+t75+t76+t77+2*f[2]+f[3]*f[3]%P)%P;

  /* Compute u_(-D1+D2) = x^3 + c1 x^2 + c2 x + c3 */
  int t82 = u1[2]*u2[2]%P;
  int t83 = u1[2]*u2[1]%P;
  int t84 = u1[1]*u2[2]%P;
  int t85 = (u1[1] + u2[1]+u1[0]+u2[0]+t82+t83+t84%P)*(1+t79+3*t59-u1[2]-u2[2]%P)%P;
  int t86 = (u1[0]+u2[0]+t83+t84%P)*(t79+3*t59-u1[2]-u2[2]%P)%P;
  c1 = t79 + 3*t59-u1[2]-u2[2]%P;
  int t87 = c1*(u1[2] + u2[2])%P;
  c2 = t80 + 3*t60 + 3*s3 - u1[1]-u2[1]-t82-t87 %P;
  int t88 = c2*(u1[2]+u2[2])%P;
  c3 = u1[1] + u2[1] + t68 + t81 + t82 + t86 + 3*t61 - t88 - t85%P;
  
  /* Compute res(t-s^2, u_(-D1-D2), x) */
  int t89 = c3*t64%P;
  int t90 = c1*t64%P;
  int t91 = c2*t64%P;
  int t92 = c2*(t65-s4)%P;
  int t93 = c1*(t66+s4+s6-s5)%P;
  int t94 = c3*(t66 + s4 +s6-s5)%P;
  int t95 = c2*(t67-s6)%P;
  int t96 = c3*(t65-s4)%P;
  int t97 = c1*(t67-s6)%P;
  int s8 = (t89 + s6-t67%P)*(t89 + s6-t67%P)%P;
  int s9 = (t91+s5-t66-s4-s6%P)*(t91+s5-t66-s4-s6%P)%P;
  int t98 = (t94 - t95)*(t90+s4-t65)%P;
  int t99 = (s8-t98)*(t89+t92+s6-t67-t93)%P;
  int t100 = (t96-t97)*(t90-t65+s4)%P;
  int t101 = (t91+s5-t66-s4-s6)*(t89+s6-t67)%P;
  int t102 = (t96-t97)*(t100-2*t101)%P;
  int t103 = s9*(t94-t95)%P;
  int res2 = t99+t102+t103%P;
  int t104 = (t90+s4-t65)*(t92+t89+s6-t93-t67)%P;
  int j0 = t104-s9%P;
  int t105 = c1*j0%P;
  int t106 = c1*(t100-t101)%P;
  int t107 = c2*j0%P;
  int t108 = c3*(t66+s4+s6-s5)%P;
  int t109 = (t108 - t95)*(t90+s4-t65)%P;
  int j1 = t105+t101-t100%P;
  int j2 = t107+t109-t106-s8%P;
  int t110 = t62*(t65+t66)%P;
  int t111 = (t62+t63)*(t66+t67)%P;
  int t112 = t63*(t65+t67)%P;
  int t113 = t63*t67%P;
  int t114 = (t62+t63)*(t66+t67)%P;
  int t115 = c1*(1-t69)%P;
  int t116 = c1*(t115+t71+t110-f[3]+t111-t69-t115-t71-t110%P)%P;
  int t117 = c2*(1-t69)%P;
  int t118 = (c2+c3)*(1+f[3]+t111-t69-t115-t71-t110%P)%P;
  int t119 = c3*(t115+t71+t110-f[3]-t111%P)%P;
  int t120 = j0*(t116 + f[2] + t113-t117-t112-t111%P)%P;
  int t121 = (j0+j1)*(t116+f[2]+f[1]+2*t113-t112-t114-t118-t119%P)%P;
  int t122 = j1*(f[1]+t111+t113+t117-t114-t118-t119%P)%P;
  int t123 = (j0+j2)*(t116+f[2]+f[0]+t119-t112-t117-t111%P)%P;
  int t124 = j2*(f[0]+t119 - t113)%P;
  int t125 = (j1+j2)*(f[1]+f[0] + t111 + t117 + t119-t114-t118-t119%P)%P;
  int t126 = c1*t120%P;
  int t127 = c2*t120%P;
  int t128 = c1*(t126 + t120 + t122 - t121%P)%P;
  int t129 = (c2+c3)*(t121-t126-t122)%P;
  int t130 = c3*(t126+t120+t122-t121)%P;

  /*Compute v_(D1+D2) */
  int t131 = res2*(t128+t123+t122-t127-t120-t124%P)%P;
  int i2 = inverse(t131);
  int t132 = i2*(t128 + t123+t122-t127-t120-t124%P)%P;
  int t133 = t132*(t128 + t123+ t122-t127-t120-t124%P)%P;
  int t134 = t132*(t125 + t127+ t122-t127-t129-t130%P)%P;
  int t135 = t132*(t130 + t124%P)%P;
  int vv2 = -t133%P;
  int vv1 = -t134%P;
  int vv0 = -t135%P;

  /*Compute u_(D1+D2) */
  int s10 = res2*res2%P;
  int t136 = i2*s10%P;
  int s11 = t136*t136%P;
  int t137 = t136*s11%P;
  int t138 = t136*t134%P;
  int s12 = t138*t138%P;
  int t139 = t136*t135%P;
  int t140 = t138*(s12+6*t139%P)%P;
  int t141 = t137*f[3]%P;
  int t142 = c1*(3*t138-c1)%P;
  int d1 = (3*t138 - c1)%P;
  int d2 = (3*t139 + 3*s12+t137-c2-t142) %P;
  int t143 = c1*d2%P;
  int t144 = c2*(3*t138-c1)%P;
  int d3 = (t140 + t141-c3-t143-t144) %P;

  printf("D_1 + D_2 = div(u,v) = div(x^3 + %d x^2 + %d x + %d, %d x^2 + %d x + %d\n", d1, d2, d3, vv2, vv1, vv0);

}

int main(int argc, char *argv[]) 
{

  /* The two divisors. Di = [ui, vi]*/
  int i;

  vsub[0] = (v1[0] - v2[0])%P;
  vsub[1] = (v1[1] - v2[1])%P;
  vsub[2] = (v1[2] - v2[2])%P;

  resultant();
  cubic();

}
