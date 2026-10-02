def _observe_with_deadline(context, live, label, *, final=False):
    end = float(live.get("case_deadline_epoch", time.time() + 300))
    remaining = min(end, float(context["deadline_epoch"])) - time.time()
    _require(remaining > 0, "Case or total deadline before native observer")
    old_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _alarm)
    old_timer = signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        return observe_price(context, live, label, final=final)
    finally:
        signal.setitimer(signal.ITIMER_REAL, *old_timer)
        signal.signal(signal.SIGALRM, old_handler)
