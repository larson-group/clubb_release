"""Tuning implementations and their shared BaseTuningStrategy interface.

algorithms.json lists the available classes, typed options and requirements.
Import an implementation from its module, or use tuner.algorithm_registry for
catalog lookup and construction through tuner.tuning_strategy. Importing this
package does not load strategy implementations or model runtimes.
"""
