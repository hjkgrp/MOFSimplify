# MOFSimplify 

The MOF website for property prediction and community engagement.

### How to run mofSimplify (tested on Python 3.8.5):
- Install Flask, molSimplify, and any other necessary dependencies.
- Run `app.py` to start a server instance.
- Go to `http://localhost:8000/` in your browser (or whatever address `app.py` prints)
- Quit and re-run `app.py` every time you make changes to the frontend or backend.

### Structure of mofSimplify:
- Backend: `app.py`
- Frontend: `index.html`
  - Dependencies: Contained in `libraries/` folder (including Bootswatch theme).
  - HTML: Any lines inside the `<body>` tag.
  - JavaScript: Any lines inside the `<script>` tag.
  - CSS: Any lines inside the `<style>` tag.
  - Interactive elements: Look at the Javascript code inside the `<script>` tag at the bottom of the file.
