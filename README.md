# MOFSimplify 

The MOF website for property prediction and community engagement.
Tested on Chrome Version 93.0.4577.82, Safari Version 15.0, Firefox 92.0, and Edge 94.0.992.37.

### How to run MOFSimplify locally (tested on Python 3.8.5):
- Install Flask, molSimplify, and any other necessary dependencies.
- Run `python app.py` to start a server instance.
- Go to `http://localhost:8000/` in your browser (or whatever address `app.py` prints)
- Quit and re-run `python app.py` every time you make changes to the frontend or backend.

### Structure of MOFSimplify:
- Backend: `app.py`
- Frontend: `index.html`
  - Dependencies: Contained in `libraries/` folder (including Bootswatch theme).
  - HTML: Any lines inside the `<body>` tag.
  - JavaScript: Any lines inside the `<script>` tag.
  - CSS: Any lines inside the `<style>` tag.
  - Interactive elements: Look at the Javascript code inside the `<script>` tag at the bottom of the file.
