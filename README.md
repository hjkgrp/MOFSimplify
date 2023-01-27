# MOFSimplify 

The MOF website for property prediction and community engagement. Available at [https://mofsimplify.mit.edu/](https://mofsimplify.mit.edu/)

Tested on Chrome Version 93.0.4577.82, Safari Version 15.0, Firefox 92.0, and Edge 94.0.992.37. 

### How to run MOFSimplify locally (tested on Python 3.8.5):
- Install Flask, molSimplify, and any other necessary dependencies (see [environments/environment.yml](environments/environment.yml)).
- Run `python app.py` to start a server instance.
  - Will need to comment out MongoClient code in `app.py` if you don't have permission to access the MOFSimplify databases. 
- Go to `http://localhost:8000/` in your browser (or whatever address `app.py` prints).
- Quit and re-run `python app.py` every time you make changes to the frontend or backend.

### Structure of MOFSimplify:
- Backend: `app.py`
  - Contains the code that generates stability predictions on MOFs.
  - Contains the code that gets a MOF's components.
  - Contains the code that generates MOFs from building blocks.
  - Contains the code that gets information on a MOF's latent space nearest neighbor.
  - Contains the code that sends information to MOFSimplify's databases.
- Frontend: `index.html`
  - Allows the user to request and see analysis on MOFs of their choosing.
  - Allows the user to give feedback and upload information about new MOFs.
  - Sends MOF information to the backend for analysis, and receives the analysis from the backend.
  - Contains the code that visualizes MOFs.
  
  - Dependencies: Contained in `libraries/` folder (including Bootswatch theme).
  - HTML: Any lines inside the `<body>` tag.
  - JavaScript: Any lines inside the `<script>` tag.
  - CSS: Any lines inside the `<style>` tag.
  - Interactive elements: Look at the Javascript code inside the `<script>` tag at the bottom of the file.

So, the most important files are `index.html` and `app.py`
