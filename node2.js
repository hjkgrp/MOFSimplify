// this file hosts the building block generated MOF for use in index.html
// https://coderrocketfuel.com/article/how-to-serve-static-files-using-node-js-and-express

const express = require("express")
const path = require("path")

const PORT = process.env.PORT || 5002

const app = express()

app.use("/temp_file_creation", express.static(path.join(__dirname, "temp_file_creation")))

app.listen(PORT, function () {
  console.log(`Express server listening on port ${PORT}`)
})