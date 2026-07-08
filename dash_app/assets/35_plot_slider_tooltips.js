(function () {
  window.dccFunctions = window.dccFunctions || {};

  function cleanNumber(value) {
    if (!Number.isFinite(value)) {
      return value;
    }
    return Number.parseFloat(value.toPrecision(3)).toString();
  }

  function formatSeconds(value) {
    var seconds = Number(value);
    var absSeconds = Math.abs(seconds);
    if (!Number.isFinite(seconds)) {
      return value;
    }
    if (absSeconds >= 3600) {
      return cleanNumber(seconds / 3600) + "h";
    }
    if (absSeconds >= 60) {
      return cleanNumber(seconds / 60) + "m";
    }
    return cleanNumber(seconds) + "s";
  }

  window.dccFunctions.formatPlotSeconds = formatSeconds;
  window.dccFunctions.formatPlotMinutes = function (value) {
    return formatSeconds(Number(value) * 60);
  };
}());
