"""Display-only layout corrections for the original seventeen checkpoint plots."""
import sys,os,json,hashlib
from pathlib import Path
import numpy as np
root=Path(os.environ["E5F_DIAGNOSTIC_SOURCE_ROOT"])
sys.path[:0]=[str(root/"code/model/tools"),str(root/"code/model")]
import run_e5f_joint_nested_finalize as f
from matplotlib.figure import Figure
selected=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
f.configure_policy_model()
receipt=f.adapter.read_json(selected/"case_receipt.json")
summary=f.adapter.read_json(selected/"summary.json")
assert f.calibration.code_fingerprint_contract(f.solver)["bundle_sha256"] == summary["code_fingerprints"]["bundle_sha256"]
f.adapter.verify(selected/"dated_state.pkl.gz",receipt["artifact_sha256"]["dated_state.pkl.gz"])
packet=f.audit.load_checkpoint(selected/"dated_state.pkl.gz")
original=Figure.savefig;checks={}
def data_hash(fig):
    data=[]
    for ax in fig.axes:
        data.append(dict(lines=[(np.asarray(x.get_xdata()).tolist(),np.asarray(x.get_ydata()).tolist()) for x in ax.lines],
            bars=[list(x.get_bbox().bounds) for x in ax.patches if hasattr(x,"get_bbox")],
            collections=[np.asarray(x.get_offsets()).tolist() for x in ax.collections],
            xlim=list(ax.get_xlim()),ylim=list(ax.get_ylim()),title=ax.get_title(),xlabel=ax.get_xlabel(),ylabel=ax.get_ylabel()))
    return hashlib.sha256(json.dumps(data,sort_keys=True,default=str).encode()).hexdigest()
def display_save(fig,path,*args,**kwargs):
    before=data_hash(fig)
    if not getattr(fig,"_review_layout_done",False):
        handles=[]; names=[]
        for ax in fig.axes:
            leg=ax.get_legend()
            if leg is not None:
                handles.extend(leg.legend_handles)
                names.extend(t.get_text() for t in leg.get_texts())
                leg.remove()
            labels=ax.get_xticklabels()
            if len(labels)>8 and max((len(x.get_text()) for x in labels),default=0)>6:
                for label in labels:
                    label.set_rotation(55);label.set_ha("right");label.set_fontsize(8)
        if handles:
            fig.canvas.draw()
            bottom=fig.get_tightbbox(fig.canvas.get_renderer()).y0/fig.get_size_inches()[1]
            fig.legend(handles,names,loc="upper center",bbox_to_anchor=(.5,min(0.,bottom)-.035),ncol=min(5,len(names)),fontsize=10 if len(fig.axes)>2 else 8,frameon=False)
            assert [x.get_text() for x in fig.legends[-1].get_texts()]==names
        fig._review_layout_done=True
    after=data_hash(fig);assert before==after,"Display edit changed plotted data or axes"
    kwargs["bbox_inches"]="tight"
    original(fig,path,*args,**kwargs)
    path=Path(path)
    if path.suffix==".png":
        origin=selected/"standard_diagnostics"/path.name
        f.adapter.verify(origin,receipt["artifact_sha256"]["standard_diagnostics/"+path.name])
        checks[path.name]=dict(original_sha256=f.adapter.digest(origin),display_sha256=f.adapter.digest(path),data_hash=after,data_unchanged=True)
Figure.savefig=display_save
try:f.audit.standard_diagnostics(packet,out,validate_production_young=False)
finally:Figure.savefig=original
assert len(checks)==17
f.adapter.write_json(out/"display_manifest.json",dict(status="same_plot_data_legibility_edits",selected_summary_sha256=f.adapter.digest(selected/"summary.json"),graphs=checks,changes="Crowded legends moved below axes; crowded numeric x labels rotated. Same seventeen plots, axis limits, labels and numerical artist data.",production_promoted=False))
