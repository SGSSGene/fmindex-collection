document$.subscribe(function() {
  var headings = document.querySelectorAll("th")
  headings.forEach(function(th) {
    th.setAttribute("data-sort-method", 'number')
  })

  var tables = document.querySelectorAll("article table:not([class])")
  tables.forEach(function(table) {
    new Tablesort(table)
  })
})
