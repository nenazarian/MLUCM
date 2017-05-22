#Major modifications to the MLUCM
in subroutine 'drag_length', parameters 'bx,by,wx,wy' are passed on so that they are not 'hard-coded' anymore
parameterization of Lk:
Changed parameterization from Leps to Lk
parameterization updated based on LkM that includes dispersive stress
The ck factors taken out and instead Ck.LkM (and similarly Leps/Ceps) is parameterized
2nd degree polynomial parameterization added for Ck.LkM above building height, a1,a2, a3 parameters introduced and calculated based on Lp
Ck.LkM below building height is calculated based on Lp (4th order polynomial!)
Parameterization of Drag:
CDeq (same in x and y) updated based on Lp (and not Lf_x and Lf_y)
by=bx assumption is taken out (how important is that if bx by is correctly passed on)
