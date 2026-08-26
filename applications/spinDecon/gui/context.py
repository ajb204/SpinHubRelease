"""GUI access to application-scoped dependencies during legacy migration.

New workspaces should receive :class:`ApplicationContext` explicitly.  Older
Frames still receive ``tabOne``; ``context_for`` provides a single temporary
bridge so those Frames no longer need to infer project/data ownership through
parent chains.
"""
from __future__ import annotations


def context_for(*objects):
    """Return the first ApplicationContext reachable from legacy GUI objects."""
    seen = set()
    queue = [obj for obj in objects if obj is not None]
    while queue:
        obj = queue.pop(0)
        ident = id(obj)
        if ident in seen:
            continue
        seen.add(ident)
        context = getattr(obj, "app_context", None)
        if context is not None:
            return context
        for name in ("parent", "Parent"):
            candidate = getattr(obj, name, None)
            if candidate is not None and candidate is not obj:
                queue.append(candidate)
        get_parent = getattr(obj, "GetParent", None)
        if callable(get_parent):
            try:
                candidate = get_parent()
            except Exception:
                candidate = None
            if candidate is not None and candidate is not obj:
                queue.append(candidate)
    return None


def project_for(*objects):
    context = context_for(*objects)
    if context is not None and context.project is not None:
        return context.project
    for obj in objects:
        state = getattr(obj, "state", None)
        if state is not None:
            return state
    return None


def data_for(*objects):
    context = context_for(*objects)
    if context is not None and context.data is not None:
        return context.data
    for obj in objects:
        store = getattr(obj, "store", None) or getattr(obj, "data_store", None)
        if store is not None:
            return store
    return None
