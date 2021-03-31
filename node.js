// this file hosts the example mof cif file for use in Load example MOF in index.html
// https://coderrocketfuel.com/article/how-to-serve-static-files-using-node-js-and-express

const express = require("express")
const path = require("path")

const PORT = process.env.PORT || 5001

const app = express()

app.use("/mof_examples", express.static(path.join(__dirname, "mof_examples")))

app.listen(PORT, function () {
  console.log(`Express server listening on port ${PORT}`)
})