# Target architecture

Dependency direction:

    app -> gui/workflow/project/domain/analysis/processing
    gui -> workflow/project/domain/analysis
    workflow -> project/domain/analysis/processing
    project -> domain
    analysis -> domain
    processing -> domain
    domain -> standard/third-party scientific libraries only

Scientific/domain modules must not depend on wx or GUI modules.

`Frames/` is frozen: new GUI code belongs under `gui/`. Existing Frames modules are migrated incrementally only after their dependencies are made explicit.
