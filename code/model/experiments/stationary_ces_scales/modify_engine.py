from pathlib import Path
root=Path(__file__).parent/'refactor_lab/engine'
p=root/'kernels.py';s=p.read_text()
helpers='''
@njit(cache=True)
def ces_flow(c, h, alpha, oms, eta, ec, eh):
    if c <= 1e-10 or h <= 1e-10:
        return -1e10
    rho = (eta - 1.0) / eta
    if abs(rho) < 1e-8:
        X = (c/ec)**alpha * (h/eh)**(1.0-alpha)
    else:
        X = (alpha*(c/ec)**rho + (1.0-alpha)*(h/eh)**rho)**(1.0/rho)
    return X**oms/oms


@njit(cache=True)
def ces_marginal_c(c, h, alpha, oms, eta, ec, eh):
    if c <= 1e-10:
        return 1e300
    rho = (eta-1.0)/eta
    if abs(rho)<1e-8:
        X = (c/ec)**alpha*(h/eh)**(1.0-alpha)
        return alpha*X**oms/c
    Q = alpha*(c/ec)**rho+(1.0-alpha)*(h/eh)**rho
    return alpha/ec*(c/ec)**(rho-1.0)*Q**((oms-rho)/rho)


@njit(cache=True)
def ces_renter_ratio(rent, alpha, eta, ec, eh):
    rho = (eta-1.0)/eta
    return (((1.0-alpha)/alpha)*(ec/eh)**rho/rent)**eta


@njit(cache=True)
def ces_renter_allocation(S, rent, hmax, alpha, oms, eta, ec, eh):
    if S <= 1e-10:
        return -1e10, 0.0, 0.0
    ratio = ces_renter_ratio(rent, alpha, eta, ec, eh)
    c = S/(1.0+rent*ratio)
    h = ratio*c
    if h > hmax:
        h = hmax
        c = S-rent*h
    return ces_flow(c,h,alpha,oms,eta,ec,eh),c,h


@njit(cache=True)
def ces_saving_mu(bp, resources, rent, hmax, alpha, oms, eta, ec, eh, owner_cost, owner_service, owner):
    S=resources-bp
    if owner:
        return ces_marginal_c(S-owner_cost,owner_service,alpha,oms,eta,ec,eh)
    _,c,h=ces_renter_allocation(S,rent,hmax,alpha,oms,eta,ec,eh)
    return ces_marginal_c(c,h,alpha,oms,eta,ec,eh)

'''
s=s.replace('@njit(cache=True)\ndef interp_scalar',helpers+'\n@njit(cache=True)\ndef interp_scalar',1)
# All existing objectives gain optional CES arguments; CD branch remains unchanged.
s=s.replace('hbc=0.0, w0=0.0, w1=0.0, hk=6.0, wedge_on=0):','hbc=0.0, w0=0.0, w1=0.0, hk=6.0, wedge_on=0, ces_eta=0.0, ces_ec=1.0, ces_eh=1.0):')
needle='    if wedge_on != 0:\n        u_flow, _, _ = renter_wedge_flow'
s=s.replace(needle,'''    if ces_eta > 0.0:
        u, _, _ = ces_renter_allocation(Rv-bp,ri,hRmax,alpha,oms,ces_eta,ces_ec,ces_eh)
        if u <= -1e9:
            return -1e10
        return u+pc+beta*interp_scalar(bg,Vbar,bp)
'''+needle,1)
s=s.replace('def eval_owner_scalar(bp, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, es):','def eval_owner_scalar(bp, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, es, ces_eta=0.0, ces_ec=1.0, ces_eh=1.0, ces_service=1.0):')
s=s.replace('    ct = Rv - oc - cb_c - bp\n    if ct <= 1e-10:', '''    ct = Rv - oc - cb_c - bp
    if ces_eta > 0.0:
        u=ces_flow(ct,ces_service,alpha,oms,ces_eta,ces_ec,ces_eh)
        if u <= -1e9: return -1e10
        return u+pc+beta*interp_scalar(bg,Vbar,bp)
    if ct <= 1e-10:''',1)
s=s.replace('def _renter_value(bp, r, Rv, Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, es):','def _renter_value(bp, r, Rv, Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, es, ces_eta=0.0, ces_ec=1.0, ces_eh=1.0):')
s=s.replace('    # eval_renter_scalar, wedge_on == 0 branch, same operation order.','''    if ces_eta > 0.0:
        u,_,_=ces_renter_allocation(Rv-bp,ri,hRmax,alpha,oms,ces_eta,ces_ec,ces_eh)
        if u <= -1e9: return -1e10
        return u+pc+beta*_interp_ranked(bg,Vbar,bp,r)
    # eval_renter_scalar, wedge_on == 0 branch, same operation order.''',1)
s=s.replace('def _owner_value(bp, r, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, es):','def _owner_value(bp, r, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, es, ces_eta=0.0, ces_ec=1.0, ces_eh=1.0, ces_service=1.0):')
start=s.index('def _owner_value');end=s.index('\n\n@njit',start)
part=s[start:end];part=part.replace('    if ct <= 1e-10:','''    if ces_eta > 0.0:
        u=ces_flow(ct,ces_service,alpha,oms,ces_eta,ces_ec,ces_eh)
        if u <= -1e9: return -1e10
        return u+pc+beta*_interp_ranked(bg,Vbar,bp,r)
    if ct <= 1e-10:''',1);s=s[:start]+part+s[end:]
s=s.replace('hmax, alpha, oms, beta, es, owner_cost, owner_K, owner):','hmax, alpha, oms, beta, es, owner_cost, owner_K, owner, ces_eta=0.0, ces_ec=1.0, ces_eh=1.0, ces_service=1.0):',1)
start=s.index('def exhaustive_saving_scalar');end=s.index('\n\n@njit',start);part=s[start:end]
part=part.replace('    nb = bg.size','''    if ces_eta > 0.0:
        if hb != 0.0 or cb != 0.0 or es != 1.0:
            raise ValueError("CES forbids floors and utility multipliers")
        ratio=ces_renter_ratio(rent,alpha,ces_eta,ces_ec,ces_eh)
        cap=hmax*(1.0+rent*ratio)/ratio
    nb = bg.size''')
part=part.replace('owner_K,alpha,oms,beta,es)','owner_K,alpha,oms,beta,es,ces_eta,ces_ec,ces_eh,ces_service)')
part=part.replace('hcap,Kr,alpha,oms,beta,es)','hcap,Kr,alpha,oms,beta,es,ces_eta,ces_ec,ces_eh)')
part=part.replace('        if owner:\n            optimal_c=', '''        if ces_eta > 0.0:
            # Within each linear-continuation interval the material utility
            # is concave. Its marginal utility rises with saving; a bracketed
            # FOC is the unique interior maximum, otherwise endpoints suffice.
            left=x;right=points[k+1];target=beta*slope
            ml=ces_saving_mu(left,resources,rent,hmax,alpha,oms,ces_eta,ces_ec,ces_eh,owner_cost,ces_service,owner)
            mr=ces_saving_mu(right,resources,rent,hmax,alpha,oms,ces_eta,ces_ec,ces_eh,owner_cost,ces_service,owner)
            if not (ml < target < mr):
                continue
            for iteration in range(44):
                midroot=(left+right)/2.0
                mu=ces_saving_mu(midroot,resources,rent,hmax,alpha,oms,ces_eta,ces_ec,ces_eh,owner_cost,ces_service,owner)
                if mu < target: left=midroot
                else: right=midroot
            candidate=(left+right)/2.0
        elif owner:
            optimal_c=''',1)
s=s[:start]+part+s[end:]
# Full block optional args, per-cell scale, propagating every scalar call.
s=s.replace('    fixed_renter_floor=-np.inf,\n):','    fixed_renter_floor=-np.inf,\n    ces_eta=0.0, ces_ec_v=None, ces_eh_v=None,\n):',1)
s=s.replace('    due_death_floor=-np.inf,\n):','    due_death_floor=-np.inf,\n    ces_eta=0.0, ces_ec_v=None, ces_eh_v=None,\n):',1)
for fname in ('full_renter_block_kernel','full_owner_block_kernel'):
 start=s.index('def '+fname);end=s.index('\n\n@njit',start);part=s[start:end]
 part=part.replace('        es = esc_v[c]','''        es = esc_v[c]
        ec = ces_ec_v[c] if ces_ec_v is not None else 1.0
        eh = ces_eh_v[c] if ces_eh_v is not None else 1.0
        if ces_eta > 0.0 and (cbc != 0.0 or hbc != 0.0 or es != 1.0):
            raise ValueError("CES full block received a floor or material multiplier")''',1)
 if fname=='full_renter_block_kernel':
  part=part.replace('0.0, 0.0, False)','0.0, 0.0, False, ces_eta, ec, eh)')
  part=part.replace('hk, wedge_on)','hk, wedge_on, ces_eta, ec, eh)')
  part=part.replace('            if wedge_on != 0:\n                uw,', '''            if ces_eta > 0.0:
                _,ct,ht=ces_renter_allocation(Rvb-bp_best,ri,hR_max,al,oms,ces_eta,ec,eh)
                co[b,c]=ct
                ho[b,c]=ht
                continue
            if wedge_on != 0:
                uw,''',1)
 else:
  part=part.replace('oc, Ko_c, True)','oc, Ko_c, True, ces_eta, ec, eh, ht_c)')
  part=part.replace('Ko_c, al, oms, beta, es)','Ko_c, al, oms, beta, es, ces_eta, ec, eh, ht_c)')
  part=part.replace('            if exact_allocation_output and exhaustive_saving and ct > 1e-10:','            if (ces_eta > 0.0 or (exact_allocation_output and exhaustive_saving)) and ct > 1e-10:')
 s=s[:start]+part+s[end:]
p.write_text(s)
# Shared CES scales use literal children at home and retain direct-child benefit.
p=root/'shared.py';s=p.read_text();s=s.replace('    apply_child_preferences(P, alpha_bar, psi_v, escale)','''    apply_child_preferences(P, alpha_bar, psi_v, escale)
    ces_ec=np.ones_like(escale)
    ces_eh=np.ones_like(escale)
    if bool(getattr(P,"ces_enabled",False)):
        if not independent_child_maturation_active(P):
            raise ValueError("CES requires independent current-child counts")
        if (np.any(c_bar != 0.0) or np.any(h_bar != 0.0)
                or bool(getattr(P,"compensated_child_housing_shares",False))
                or float(P.delta_alpha)!=0.0 or float(P.delta_alpha_jump)!=0.0
                or float(P.sigma)!=2.0 or float(P.alpha_cons)!=.733):
            raise ValueError("CES primitive contract drift: no floors/A/share variation, sigma2, alpha.733")
        escale[:]=1.0
        alpha_bar[:]=P.alpha_cons
        for nn in range(P.n_parity):
            for cs in range(P.n_child_states):
                m=cs if cs<=nn else 0
                ces_ec[nn,cs]=((2.0+.7*m)/2.0)**.7
                ces_eh[nn,cs]=ces_ec[nn,cs]*(1.0+float(P.lambda_housing)*(m>0))''',1)
s=s.replace('        escale_flat=escale.reshape(1, nc, order="F"),','        escale_flat=escale.reshape(1, nc, order="F"),\n        ces_ec_flat=ces_ec.reshape(1,nc,order="F"),\n        ces_eh_flat=ces_eh.reshape(1,nc,order="F"),',1);p.write_text(s)
p=root/'household.py';s=p.read_text();s=s.replace('        esc_v=np.ascontiguousarray(SD.escale_flat.reshape(-1)),','        esc_v=np.ascontiguousarray(SD.escale_flat.reshape(-1)),\n        ces_eta=float(getattr(P,"ces_eta",0.0)) if bool(getattr(P,"ces_enabled",False)) else 0.0,\n        ces_ec_v=np.ascontiguousarray(SD.ces_ec_flat.reshape(-1)),\n        ces_eh_v=np.ascontiguousarray(SD.ces_eh_flat.reshape(-1)),',1)
s=s.replace('                float(renter_floor[0]) if fixed_credit else -np.inf,','                float(renter_floor[0]) if fixed_credit else -np.inf,\n                ctx.ces_eta, ctx.ces_ec_v, ctx.ces_eh_v,',1)
s=s.replace('                if natural_credit:\n                    Vo_nc[:, natural_dead] = -1e10\n            else:', '                if natural_credit:\n                    Vo_nc[:, natural_dead] = -1e10\n            else:',1)
s=s.replace('native_due_death_floor(P, j, ctx.current_prices[i], P.H_own[ten - 1]) if due_stay else -np.inf,','native_due_death_floor(P, j, ctx.current_prices[i], P.H_own[ten - 1]) if due_stay else -np.inf,\n                    ctx.ces_eta, ctx.ces_ec_v, ctx.ces_eh_v,',1)
s=s.replace('    use_full_kernel = NUMBA_AVAILABLE and bool(getattr(P, "use_full_kernel", True)) and interp_method == "linear"','''    use_full_kernel = NUMBA_AVAILABLE and bool(getattr(P, "use_full_kernel", True)) and interp_method == "linear"
    if bool(getattr(P,"ces_enabled",False)) and not use_full_kernel:
        raise ValueError("CES requires the integrated full indexed kernel; CD fallback forbidden")''',1)
p.write_text(s)
