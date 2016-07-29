# remove NOTE about no visible binding for global variable during R CMD check --
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("value", "Pop", "variable", ".", "newid", "ID", "PopNames", "Name", "Names", "Group", "major")
  )
}
